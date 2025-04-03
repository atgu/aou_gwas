source('~/CHARR/R/constants.R')
library(tidyverse)
library(RColorBrewer)
# install.packages('~/Downloads/maptools_1.1-8.tar.gz', repos=NULL, type='source')
# .tar.gz file downloaded from https://cran.r-project.org/src/contrib/Archive/maptools/
library(maptools)
library(cowplot)
library(ggpubr)

setwd('~/Dropbox (Partners HealthCare)/aou/pca/data/')

pop_colors = c('AFR' = '#941494', 'AMR' = '#ED1E24', 'SAS' = '#FF9912', 'EAS' = '#108C44', 'EUR' = '#6AA5CD', 'MID' = '#EEA9B8', 'OCE' = '#a6761d', 'UKBB' = 'black', 'Other' = '#ABB9B9')
pops <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')
tag = '_pruned'

meta_info <- read_delim('aou_meta_info_250k.txt.bgz', delim='\t')
#### Load pop-specific PCAs
aou_pc_full <- read_delim(paste0('aou_sas_pca', tag,'_full_scores_250k.txt.bgz'), delim='\t') %>%
  filter(!is.na(PC1)) %>%
  mutate(pop = 'SAS',
         s = as.character(s))
for(pop1 in pops[1:5]){
  # if(pop1 == 'eur') next
  aou_pc <- read_delim(paste0('aou_', pop1,'_pca',tag,'_full_scores_250k.txt.bgz'), delim='\t')%>%
    filter(!is.na(PC1)) %>%
    mutate(pop = toupper(pop1),
           s = as.character(s))
  aou_pc_full <- rbind(aou_pc_full, aou_pc)
}

aou_pc_full <- aou_pc_full %>%
  merge(., meta_info %>% select(person_id, global_PC1 = PC1, race = race), by.x = 's', by.y='person_id')

## Global PCA
plot_global_pca <- function(pcs, pop_color, pop_shape=NA, first_pc='PC1', second_pc='PC2', legend_name='pop') {
  if(is.na(pop_shape)) {
    pca_pop <- ggplot(pcs, aes_string(x=first_pc, y=second_pc, color=legend_name)) + geom_point(alpha=0.75)
  } else {
    pca_pop <- ggplot(pcs, aes_string(x=first_pc, y=second_pc, color=legend_name, shape=legend_name)) + geom_point(alpha=0.75)
  }
  pca_pop <- pca_pop +
    scale_color_manual(values=pop_color, name=legend_name) +
    scale_shape_manual(values=pop_shape, name=legend_name) +
    theme_classic() +
    theme(text = element_text(size=8),
          axis.text = element_text(color='black'),
          legend.text = element_text(size=10)) +
    guides(colour = guide_legend(nrow = 1)) +
    themes


  x_lim = ggplot_build(pca_pop)$layout$panel_scales_x[[1]]$range$range
  y_lim = ggplot_build(pca_pop)$layout$panel_scales_y[[1]]$range$range

  return(list(pca_pop, x_lim, y_lim))
}

for(pop1 in pops){
  global_pcs <- list()
  # if( pop1 != 'eur') next
  # if( ! pop1 %in% c('amr', 'afr')) next
  for(i in 1:25){
    global_pcs[[i]] <- plot_global_pca(aou_pc_full %>% filter(pop == toupper(pop1)), pop_colors, first_pc = paste0('PC', 2*i-1), second_pc = paste0('PC', 2*i), legend_name = 'pop')[[1]]
    # png(paste0('~/Desktop/aou_pca_', 2*i-1, '_', 2*i, '_old_pruned.png'), width=4, height=3, units = 'in', res = 300)
    # print(global_pcs[[i]])
    # dev.off()
  }
  global_pcs_sub = ggpubr::ggarrange(plotlist = global_pcs[1:10],
                                   ncol = 5, nrow=2, common.legend = TRUE, labels=LETTERS[1:10] )
  png(paste0('~/Desktop/aou_pca_sub_', pop1, tag,'_by_race.png'), width=9, height=4, units = 'in', res = 300)
  print(global_pcs_sub)
  dev.off()

  global_pc_plot = ggpubr::ggarrange(plotlist = global_pcs,
                                   ncol = 5, nrow=5, common.legend = TRUE, labels=LETTERS[1:25] )
  png(paste0('~/Desktop/aou_pca_full_', pop1, tag,'_by_race.png'), width=9, height=9, units = 'in', res = 300)
  print(global_pc_plot)
  dev.off()
}

p <- aou_pc_full %>%
  filter(pop == 'AFR') %>%
  ggplot + aes(x = PC1, color=race) +
  labs(x = 'PC1 - AFR', y = 'Density')+
  geom_density() + themes +
  theme(legend.position = 'right') +
    guides(colour = guide_legend(ncol = 1))
png(paste0('~/Desktop/aou_pca_afr_pc1_by_race.png'), width=6, height=4, units = 'in', res = 300)
print(p)
dev.off()

## panels by pop
for(pc_id in 1:10){
  global_pcs_1 <- list()
  for(i in 1:length(pops)){
    pop1 <- toupper(pops[i])
    global_pcs_1[[i]] <- plot_global_pca(aou_pc_full %>% filter(pop == pop1), pop_colors, first_pc = paste0('PC', 2*pc_id-1), second_pc = paste0('PC', 2*pc_id))[[1]]
  }

  global_pcs_pop = ggpubr::ggarrange(plotlist = global_pcs_1,
                                     ncol = 3, nrow=2, labels=NULL )
  png(paste0('~/Desktop/aou_pca_pop_pc', 2*pc_id-1, '_', 2*pc_id, '.png'), width=6, height=4, units = 'in', res = 300)
  print(global_pcs_pop)
  dev.off()
}


aou_pc_full %>%
  filter(PC1 >-100 & PC1 <100) %>%
  dplyr::group_by(pop) %>%
  dplyr::summarize(count = n())


## Self-reported Ancestry
blues <- brewer.pal(5, 'Blues')[2:5] #1-1003
reds <- brewer.pal(6, 'Reds')[2:6] #2-2004
oranges <- brewer.pal(6, 'Oranges')[2:6] #3-3004
purples <- brewer.pal(5, 'Purples')[2:5] #4-4003
greys <- brewer.pal(4, 'Greys')[2:4] #-3, -1, 6
ethnicity_colors = c(greys[1:4], blues[1], reds[1], oranges[1], purples[1],'#33CC33', '#108C44')
names(ethnicity_colors) <- c('I prefer not to answer', 'None Indicated', 'None of these', 'Skip',
                             'White',
                             'More than one population',
                             'Asian',
                             'Black or African American',
                             'Middle Eastern or North African',
                             'Native Hawaiian or Other Pacific Islander')

ethnicity_pcs_1_2 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC1', second_pc = 'PC2', legend_name='race')
ethnicity_pcs_3_4 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC3', second_pc = 'PC4', legend_name='race')
ethnicity_pcs_5_6 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC5', second_pc = 'PC6', legend_name='race')
ethnicity_pcs_7_8 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC7', second_pc = 'PC8', legend_name='race')
ethnicity_pcs_9_10 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC9', second_pc = 'PC10', legend_name='race')
ethnicity_pcs_11_12 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC11', second_pc = 'PC12', legend_name='race')
ethnicity_pcs_13_14 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC13', second_pc = 'PC14', legend_name='race')
ethnicity_pcs_15_16 <- plot_global_pca(aou_pop, ethnicity_colors, first_pc = 'PC15', second_pc = 'PC16', legend_name='race')

ethnicity_pcs <- plot_grid(ethnicity_pcs_1_2[[1]], ethnicity_pcs_3_4[[1]], ethnicity_pcs_5_6[[1]],ethnicity_pcs_7_8[[1]],
                        ethnicity_pcs_9_10[[1]], ethnicity_pcs_11_12[[1]], ethnicity_pcs_13_14[[1]], ethnicity_pcs_15_16[[1]],
                        labels=LETTERS[1:8], nrow=2)

ethnicity_pcs = ggpubr::ggarrange(ethnicity_pcs_1_2[[1]], ethnicity_pcs_3_4[[1]], ethnicity_pcs_5_6[[1]],ethnicity_pcs_7_8[[1]],
                               ethnicity_pcs_9_10[[1]], ethnicity_pcs_11_12[[1]], ethnicity_pcs_13_14[[1]], ethnicity_pcs_15_16[[1]],
                               ncol = 4, nrow=2, common.legend = TRUE, labels=LETTERS[1:8] )
png('aou_ethnicity_pca.png', width=9, height=4, units = 'in', res = 300)
print(ethnicity_pcs)
dev.off()

ethnicity_pcs_sub = ggpubr::ggarrange(ethnicity_pcs_1_2[[1]], ethnicity_pcs_3_4[[1]], ethnicity_pcs_5_6[[1]],ethnicity_pcs_7_8[[1]],
                                     ncol = 2, nrow=2, common.legend = TRUE, labels=LETTERS[1:4] )
png('aou_ethnicity_pca_main.png', width=6, height=4, units = 'in', res = 300)
print(ethnicity_pcs_sub)
dev.off()

## PCA centroid
plot_pca_density <- function(dataset, first_pc, second_pc) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans = "log", breaks=c(1,20,400,8000,163000), name='Count',
                         colours = rev(brewer.pal(5,'Spectral'))) +
    theme_classic() +
    theme(text = element_text(size=16))
  return(pc_biplot)
}

pop_ellipse <- function(df, num_ellipses) {
  # get mean and SD of each PC among each pop
  pc_nams <- paste("PC",1:10,sep="")
  mean_pcs <- colMeans(df[,pc_nams])
  sd_pcs <- apply(df[,pc_nams],2,sd)
  # compute centroid distance for each individual
  centroid_dist <- rep(0,nrow(df))
  for(i in 1:num_ellipses) {
    centroid_dist <- centroid_dist + (df[,pc_nams[i]]-mean_pcs[i])^2/(sd_pcs[i]^2)
  }
  pop_dist <- df %>%
    mutate(centroid_dist=centroid_dist)
  return(pop_dist)
}

pop_centroid <- function(ind_dist, cutpoint0, cutpoint1) {
  pop_cut <- subset(ind_dist, centroid_dist < cutpoint0)
  p_centroid <- ggplot(pop_cut, aes(x=centroid_dist)) +
    geom_histogram(bins=50) +
    labs(title=paste0('Sample size: ', nrow(subset(ind_dist, centroid_dist < cutpoint0)), ' -> ', nrow(subset(ind_dist, centroid_dist < cutpoint1)))) +
    xlab('Centroid distance') +
    ylab('Count') +
    geom_vline(xintercept=cutpoint1) +
    theme_bw() +
    theme(text = element_text(size=16))
  return(list(p=p_centroid, pop_cut=pop_cut))
}

save_filt_plots <- function(pop_name, pop_dist, cutpoint0, cutpoint1) {
  p_centroid0 = pop_centroid(pop_dist, cutpoint0, cutpoint1)
  ggsave(paste0('../figures/',pop_name, '/aou_', pop_name, '_within_pop_centroid_nofilt.pdf'), p_centroid0$p, height=7, width=7, limitsize = FALSE)
  p2 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC1', 'PC2')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc1_2.png'), p2, limitsize = FALSE)
  p3 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC3', 'PC4')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc3_4.png'), p3, limitsize = FALSE)
  p4 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC5', 'PC6')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc5_6.png'), p4, limitsize = FALSE)
  p_centroid1 = pop_centroid(pop_dist, cutpoint1, cutpoint1)
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_centroid_filt.pdf'), p_centroid1$p, height=7, width=7, limitsize = FALSE)
  p6 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC1', 'PC2')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_filt_pc1_2.png'), p6, limitsize = FALSE)
  p7 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC3', 'PC4')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_filt_pc3_4.png'), p7, limitsize = FALSE)
  p8 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC5', 'PC6')
  ggsave(paste0('../figures/',pop_name, '/aou_',pop_name, '_within_pop_filt_pc5_6.png'), p8, limitsize = FALSE)
  my_plot=plot_grid(p_centroid0$p, p2, p3, p4, p_centroid1$p, p6, p7, p8, nrow=2)
  save_plot(paste0('../figures/aou_',pop_name, '_within_pop.png'), my_plot, base_height=10, base_width = 18, limitsize = FALSE)
  return(p_centroid1$pop_cut)
}

afr_cut <- save_filt_plots('afr', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'AFR')), 3), 100, 8)
amr_cut <- save_filt_plots('amr', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'AMR')), 3), 100, 8)
eas_cut <- save_filt_plots('eas', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'EAS')), 3), 100, 6.5)
eur_cut <- save_filt_plots('eur', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'EUR')), 5), 100, 10)
mid_cut <- save_filt_plots('mid', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'MID')), 5), 100, 5)
sas_cut <- save_filt_plots('sas', pop_ellipse(as.data.frame(aou_pc_full %>% filter(pop == 'SAS')), 3), 100, 5)

pop_cuts <- sas_cut %>%
  # bind_rows(eas_cut, mid_cut) %>%
  bind_rows(afr_cut, eas_cut, amr_cut, mid_cut, eur_cut) %>%
  select(s)

write.table(pop_cuts, 'aou_pops_centroid_pruned.tsv', row.names=F, sep='\t', quote=F)


sas_cut %>%
  # bind_rows(eas_cut, mid_cut) %>%
  bind_rows(afr_cut, eas_cut, amr_cut, mid_cut, eur_cut) %>%
  dplyr::group_by(pop) %>%
  dplyr::summarize(n = cnt())