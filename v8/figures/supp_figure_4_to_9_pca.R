library(tidyverse)
library(RColorBrewer)
install.packages('maptools_1.1-8.tar.gz', repos=NULL, type='source')
# .tar.gz file downloaded from https://cran.r-project.org/src/contrib/Archive/maptools/
library(maptools)
library(cowplot)
library(ggpubr)
packages = c('dplyr', 'optparse', 'tidyverse', 'Matrix', 'sparseMVN', 'spdep', 'qlcMatrix', 'ggplot2', 'stringr', 'ggpubr')
for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
}

#### Parameters
pop_colors = c('AFR' = '#941494', 'AMR' = '#ED1E24', 'SAS' = '#FF9912', 
               'EAS' = '#108C44', 'EUR' = '#6AA5CD', 'MID' = '#EEA9B8', 
               'OCE' = '#a6761d', 'UKBB' = 'black', 'Other' = '#ABB9B9')
pops <- c('afr', 'amr', 'eas', 'eur', 'mid', 'sas')
# pca_figure_root <- '~/Dropbox (Partners HealthCare)/aou_v8/figures/pca/'
pca_figure_root <- '~/Dropbox (Partners HealthCare)/aou_v8/figures/supplement/'
figure_index <- c('afr' = '4', 'amr' = '5', 'eas' = '6', 'eur' = '7', 'mid' = '8', 'sas' = '9')

# Example modification within your function
process_pc_scores <- function(ancestry){
  gcs_path <- paste0('gs://aou_analysis/v8/data/utils/pca/results/pca_', toupper(ancestry),'_scores_full.txt.bgz')
  local_temp_file <- tempfile(fileext = ".txt.bgz") # Create a temporary file name

  # Construct the gsutil command
  gsutil_command <- paste("gsutil cp", gcs_path, local_temp_file)

  # Execute the command (add error handling if needed)
  print(paste("Copying from GCS:", gsutil_command))
  system(gsutil_command)

  # Check if the file was copied successfully
  if (!file.exists(local_temp_file)) {
    stop(paste("Failed to copy file from GCS:", gcs_path))
  }

  # Read the local file
  print(paste("Reading local file:", local_temp_file))
  data <- readr::read_delim(local_temp_file, delim='\t') %>%
    dplyr::filter(!is.na(PC1)) %>%
    dplyr::mutate(ancestry = toupper(ancestry),
           s = as.character(s))

  # Clean up the temporary file
  unlink(local_temp_file)

  return(data)
}

aou_pc_full <- rbind(
  process_pc_scores('AFR'),
  process_pc_scores('AMR'),
  process_pc_scores('EAS'),
  process_pc_scores('EUR'),
  process_pc_scores('MID'),
  process_pc_scores('SAS')
)

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
  ggsave(paste0(pca_figure_root, pop_name, '/aou_', pop_name, '_within_pop_centroid_nofilt.pdf'), p_centroid0$p, height=7, width=7, limitsize = FALSE)
  p2 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC1', 'PC2')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc1_2.png'), p2, limitsize = FALSE)
  p3 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC3', 'PC4')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc3_4.png'), p3, limitsize = FALSE)
  p4 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC5', 'PC6')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_nofilt_pc5_6.png'), p4, limitsize = FALSE)
  p_centroid1 = pop_centroid(pop_dist, cutpoint1, cutpoint1)
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_centroid_filt.pdf'), p_centroid1$p, height=7, width=7, limitsize = FALSE)
  p6 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC1', 'PC2')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_filt_pc1_2.png'), p6, limitsize = FALSE)
  p7 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC3', 'PC4')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_filt_pc3_4.png'), p7, limitsize = FALSE)
  p8 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC5', 'PC6')
  ggsave(paste0(pca_figure_root, pop_name, '/aou_',pop_name, '_within_pop_filt_pc5_6.png'), p8, limitsize = FALSE)
  my_plot=plot_grid(p_centroid0$p, p2, p3, p4, p_centroid1$p, p6, p7, p8, nrow=2)
  save_plot(paste0(pca_figure_root, 'supp_figure', figure_index[pop_name],'_aou_',pop_name, '_within_pop.png'), my_plot, base_height=10, base_width = 18, limitsize = FALSE)
  return(p_centroid1$pop_cut)
}

afr_cut <- save_filt_plots('afr', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'AFR')), 3), 100, 10) 
# options(scipen = 999)
write_tsv(afr_cut %>% mutate(s=as.character(s)) %>% select(s), '~/Downloads/aou_afr_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_afr_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_afr_centroid_pruned.tsv"), intern=T)
# gsutil cat gs://aou_analysis/v8/data/utils/pca/results/aou_afr_centroid_pruned.tsv \
# | awk 'BEGIN{FS=OFS="\t"}
#        NR==1 {
#          for(i=1;i<=NF;i++) if($i=="s") col=i
#          print
#          next
#        }
#        NR==2 {
#          if($col=="1e+06") $col="1000000"
#        }
#        {print}
#      ' \
# | gsutil cp - gs://aou_analysis/v8/data/utils/pca/results/aou_afr_centroid_pruned.tsv


amr_cut <- save_filt_plots('amr', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'AMR')), 3), 100, 7.5) 
write_tsv(amr_cut %>% select(s), '~/Downloads/aou_amr_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_amr_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_amr_centroid_pruned.tsv"), intern=T)

eas_cut <- save_filt_plots('eas', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'EAS')), 3), 100, 8.5)
write_tsv(eas_cut %>% select(s), '~/Downloads/aou_eas_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_eas_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_eas_centroid_pruned.tsv"), intern=T)

eur_cut <- save_filt_plots('eur', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'EUR')), 5), 100, 26)
write_tsv(eur_cut %>% select(s), '~/Downloads/aou_eur_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_eur_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_eur_centroid_pruned.tsv"), intern=T)


mid_cut <- save_filt_plots('mid', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'MID')), 5), 1000, 6)
write_tsv(mid_cut %>% select(s), '~/Downloads/aou_mid_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_mid_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_mid_centroid_pruned.tsv"), intern=T)

sas_cut <- save_filt_plots('sas', pop_ellipse(as.data.frame(aou_pc_full %>% filter(ancestry == 'SAS')), 3), 100, 6)
write_tsv(sas_cut %>% select(s), '~/Downloads/aou_sas_centroid_pruned.tsv')
system(paste0("gsutil cp ~/Downloads/aou_sas_centroid_pruned.tsv gs://aou_analysis/v8/data/utils/pca/results/aou_sas_centroid_pruned.tsv"), intern=T)


pop_cuts <- sas_cut %>%
  # bind_rows(eas_cut, mid_cut) %>%
  bind_rows(afr_cut, eas_cut, amr_cut, mid_cut, eur_cut) %>%
  select(s)


related_samples <- read_delim('relatedness_flagged_samples.txt.bgz', delim = '\t')
pop_cuts %>%
  mutate(related = s %in% related_samples$sample_id) %>%
  group_by(ancestry, related) %>%
  dplyr::summarize(
    cnt = n()
  )

sas_cut %>%
  # bind_rows(eas_cut, mid_cut) %>%
  bind_rows(afr_cut, eas_cut, amr_cut, mid_cut, eur_cut) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::summarize(n = n())

# pop_cuts <- bind_rows(eas_cut, mid_cut) 

write.table(pop_cuts, 'aou_pops_centroid_pruned.tsv', row.names=F, sep='\t', quote=F)
# system(paste0("gsutil cp /home/rstudio/aou_pops_centroid_pruned.tsv gs://fc-secure-4c290511-a246-4aa5-a220-85a456bb470e/utils/pca/results/aou_pops_centroid_pruned.tsv"), intern=T)
