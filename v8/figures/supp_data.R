source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

get_freq_intervals <- function(x, endpoints, labels = NULL, right = FALSE) {
  if (is.null(labels)) {
    labels <- paste0(
      "[", head(endpoints, -1), ", ", tail(endpoints, -1), "]"
    )
  }
  cut(x,
      breaks = endpoints,
      labels = labels,
      right = right,
      include.lowest = TRUE)
}

unique_acaf <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_ACAF_signals_qced.txt.bgz'))  %>% 
  filter(Pvalue < 5e-8) %>%
  select(locus, alleles, phenoname, description) %>%
  mutate(ancestry = 'meta') 
unique_exome <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_Exome_signals_qced.txt.bgz')) %>% 
  filter(Pvalue < 5e-8)%>%
  select(locus, alleles, phenoname, description) %>%
  mutate(ancestry = 'meta')
for(ancestry in ANCESTRIES){
  ancestry = toupper(ancestry)
  print(paste0('#########', ancestry, ':'))
  print('---------ACAF associations---------')
  tmp <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_', ancestry,'_ACAF_signals_qced.txt.bgz'))  %>% 
    filter(Pvalue < 5e-8) %>%
    mutate(AF_bin = get_freq_intervals(AF_Allele2, c(0, 0.0001, 0.001, 0.01, 0.1, 1)))
  print(paste0('N associations:', nrow(tmp)))
  print(table(tmp$AF_bin))
  print(table(tmp$category))
  tmp <- tmp %>%
    select(locus, alleles, phenoname, description) %>%
    mutate(ancestry = tolower(ancestry))
  unique_acaf <- unique_acaf %>%
    merge(., tmp,  by = c('locus', 'alleles', 'phenoname', 'description'), suffix = c('', paste0('.', tolower(ancestry))), all=T)  
  
  
  print('---------Exome associations---------')
  tmp <-  read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_', ancestry,'_Exome_signals_qced.txt.bgz')) %>% 
    filter(Pvalue < 5e-8) %>%
    mutate(AF_bin = get_freq_intervals(AF_Allele2, c(0, 0.0001, 0.001, 0.01, 0.1, 1)))
  print(paste0('N associations:', nrow(tmp)))
  print(table(tmp$AF_bin))
  print(table(tmp$category))
  
  tmp <- tmp %>%
    select(locus, alleles, phenoname, description) %>%
    mutate(ancestry = tolower(ancestry))
  unique_exome <- unique_exome %>%
    merge(., tmp, by = c('locus', 'alleles', 'phenoname', 'description'), suffix = c('', paste0('.', tolower(ancestry))), all=T)  
}
unique_acaf <- unique_acaf %>%
  mutate(n_assoc = as.numeric(!is.na(ancestry.meta)) + as.numeric(!is.na(ancestry.sas)) + as.numeric(!is.na(ancestry.mid)) + 
           as.numeric(!is.na(ancestry.eur)) + as.numeric(!is.na(ancestry.eas)) + as.numeric(!is.na(ancestry.amr)) + 
           as.numeric(!is.na(ancestry.afr)) ) %>%
  filter(n_assoc > 0)%>%
  filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11')))

print(paste0('N ACAF unique associations:', nrow(unique_acaf)))

unique_exome <- unique_exome %>%
  mutate(n_assoc = as.numeric(!is.na(ancestry.meta)) + as.numeric(!is.na(ancestry.sas)) + as.numeric(!is.na(ancestry.mid)) + 
           as.numeric(!is.na(ancestry.eur)) + as.numeric(!is.na(ancestry.eas)) + as.numeric(!is.na(ancestry.amr)) + 
           as.numeric(!is.na(ancestry.afr)) ) %>%
  filter(n_assoc > 0) %>%
  filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11')))
print(paste0('N unique associations:', nrow(unique_exome)))




for(ancestry in ANCESTRIES){
  ancestry = toupper(ancestry)
  print(paste0('#########', ancestry, ':'))
  for(p in c('Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue')){
    for(maxMAF in c('0.01', '0.001', '0.0001', 'cauchy')){
      print(paste0('----------', p, '--------', maxMAF, '----------'))
      tmp <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_', ancestry,'_gene_signals_', p,'_', maxMAF,'_raw.txt.bgz'))  
      tmp$P <- tmp[[p]]
      tmp <- tmp %>% 
        filter(P < 6.7e-7)
      if(p == 'Pvalue_Burden' & maxMAF == 0.001){
        print(table(tmp$category))
      }
      print(paste0('N associations:', nrow(tmp)))
      print(table(tmp$max_MAF))
      if(maxMAF != 'cauchy'){
        print(table(tmp$annotation))
      }
    }
  }
}




for(p in c('Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue')){
  for(maxMAF in c('0.01', '0.001', '0.0001')){
    print(paste0('----------', p, '--------', maxMAF, '----------'))
    unique <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_META_gene_signals_', p,'_', maxMAF,'_raw.txt.bgz'))  
    
    unique$P <- unique[[p]]
    unique <- unique %>% 
      filter(P < 6.7e-7)
    unique <- unique %>%
      select(gene_id, gene_symbol, annotation, max_MAF, phenoname, description, P) %>%
      mutate(ancestry = 'META')
    unique_long <- unique
    
    for(ancestry in ANCESTRIES){
      ancestry = toupper(ancestry)
      print(paste0('#########', ancestry, ':'))
      
      tmp <- read_gcs_fast(paste0('gs://aou_wlu/v8_analysis/paper_preparation/summary/aou_v8_', ancestry,'_gene_signals_', p,'_', maxMAF,'_raw.txt.bgz')) 
      if(nrow(tmp) == 0){
        unique[[paste0('ancestry.', tolower(ancestry))]] = NA
        next
      } 
      tmp$P <- tmp[[p]]
      tmp <- tmp %>% 
        filter(P < 6.7e-7) %>%
        select(gene_id, gene_symbol, annotation, max_MAF, phenoname, description, P) %>%
        mutate(ancestry = ancestry)
      unique <- unique %>%
        merge(., tmp, by = c('gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phenoname', 'description'), all = T, suffix = c('', paste0('.', tolower(ancestry))))
      unique_long = rbind(unique_long, tmp)
    }
    unique <- unique %>%
      mutate(n_assoc = as.numeric(!is.na(ancestry.meta)) + as.numeric(!is.na(ancestry.sas)) + as.numeric(!is.na(ancestry.mid)) + 
                                    as.numeric(!is.na(ancestry.eur)) + as.numeric(!is.na(ancestry.eas)) + as.numeric(!is.na(ancestry.amr)) + 
                                                                                                                     as.numeric(!is.na(ancestry.afr)) ) %>%
      filter(n_assoc > 0) %>%
      filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11')))
    print(table(unique$annotation))
    print(paste0('N unique signal:', nrow(unique)))
    
    print(paste0('N unique signal:', nrow(unique_long %>%
                                            filter(!startsWith(phenoname, 'random') & !(phenoname %in% c('PP_928.1', 'PP_932.11'))) %>% select(phenoname, gene_id) %>% distinct())))
    
  }
}
## supplementary data 3: manuscript_numbers.R 




