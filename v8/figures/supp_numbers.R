source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')

######### AoU UKB phenotype mapping #########
v8_pheno_info <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/phenotype/aou_v8_phenotype_info_qced_both.tsv')
# aou_pheno <- v8_pheno_info %>%
#   filter(trait_type == 'continuous')%>%
#   mutate(description = str_to_sentence(description)) %>%
#   select(aou_phenocode = phenoname, description = description, aou_trait_type=trait_type) %>%
#   distinct()
# ukb_pheno <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/wlu_pan_ukb_phenotype_info.txt.bgz', delim ='\t')
# ukb_pheno <- ukb_pheno %>%
#   filter(trait_type %in% c('continuous', 'biomarkers') & modifier == 'irnt') %>%
#   mutate(
#     description = if_else(trait_type == 'icd10', str_replace(description, paste0(phenocode, ' '), ''), description),
#     description = if_else(phenocode == '50', 'height', description),
#     description = if_else(phenocode == 'whr', 'WHR', description),
#     description = if_else(description == 'Body mass index (BMI)', 'BMI', description),
#     description= str_to_sentence(description)
#   ) %>%
#   select(ukb_phenocode = phenocode, description=description, ukb_coding_description = description_more) %>%
#   mutate(description = if_else(is.na(description), ukb_phenocode, description)) %>%
#   distinct()
# 
# matched_pheno <- stringdist_left_join(aou_pheno, ukb_pheno, by = 'description', max_dist=2, distance_col = 'distance')
# write_csv(matched_pheno, '~/Desktop/aou_ukb_matched_quantitative_phenotype.csv') -> manual curation -> saved to /Users/wlu/Partners HealthCare Dropbox/Wenhan Lu/aou_v8/data/aou_ukb_map
genebass_pheno <- read_delim('~/Dropbox (Partners HealthCare)/aou/aou_genebass_meta/genebass_pheno_info.txt.bgz') %>%
  filter(grepl('icd', trait_type)) %>%
  select(ukb_phenocode = phenocode, description, category) %>%
  mutate(icd10  = str_split(str_replace(description, 'Date ', ''), ' ') %>% map_chr(., 1)) %>%
  select(ukb_phenocode, icd10, description)
pan_ukb_pheno <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/wlu_pan_ukb_phenotype_info.txt.bgz') %>%
  select(phenocode, description) %>%
  distinct()
phecode_info <- read_csv('~/Dropbox (Partners HealthCare)/aou/phenotype/phecodeX_info.csv')
ukb_aou_map <- read_csv('~/Dropbox (Partners HealthCare)/aou_v8/data/aou_ukb_map/aou_ukb_matched_quantitative_phenotype_v8.csv')
icd_phecode <- read_gcs_fast('gs://aou_wlu/v8_analysis/icd_phecode_map_0.8.txt.bgz') 
phecodex_phecode <- read_gcs_fast('gs://aou_wlu/v8_analysis/phecode_phecodex_map_0.8.txt.bgz')
# Join directly: each phecodeX has 1 phecode (from Hail), each phecode has 1 ICD (from Hail)
# Only conflict: multiple phecodeX may land on the same ICD — resolve once at the end
icd_phecode_clean <- icd_phecode %>%
  mutate(phecode = as.character(phecode)) %>%
  select(icd10 = ICD10, phecode, r2.icd = r2, r.icd = r)

phecodex_phecode_clean <- phecodex_phecode %>%
  mutate(phecode = as.character(Phecode)) %>%
  select(phenoname = PhecodeX, phecode, r2.phecodex = r2, r.phecodex = r)

phecodex_map <- merge(phecodex_phecode_clean, icd_phecode_clean, by = 'phecode', all.x = T) %>%
  merge(., genebass_pheno, by = 'icd10', all.x = T) %>%
  merge(., phecode_info %>% select(phenoname = phecode, category), by = 'phenoname', all.x = T) %>%
  mutate(r2_combined = r2.phecodex * r2.icd) %>%
  group_by(ukb_phenocode) %>%
  slice_max(r2_combined, with_ties = FALSE) %>%   # if multiple phecodeX map to same UKB pheno, keep best
  ungroup() %>%
  merge(., v8_pheno_info %>% select(phenoname, description) %>% distinct(), by = 'phenoname', all.x = T) %>%
  ungroup()%>%
  filter(!is.na(description.y) & !is.na(description.x))

phecode_map <- phecodex_phecode_clean %>%
  merge(., pan_ukb_pheno, by.y = 'phenocode', by.x = 'phecode', all.x = T)%>%
  merge(., phecode_info %>% select(phenoname = phecode, category), by = 'phenoname', all.x = T) %>%
  group_by(phecode) %>%
  slice_max(r2.phecodex, with_ties = FALSE) %>%   # if multiple phecodeX map to same UKB pheno, keep best
  ungroup() %>%
  merge(., v8_pheno_info %>% select(phenoname, description) %>% distinct(), by = 'phenoname', all.x = T) %>%
  ungroup()%>%
  filter(!is.na(description.y) & !is.na(description.x))
  


phecodex_map %>%
  dplyr::summarize(cnt = count(category))

head(phecodex_map)
head(ukb_aou_map)
length(unique(phecodex_map$phenoname))
length(unique(phecodex_map$icd10))
length(unique(phecodex_map$phecode))

aou_v8_mapped_ukb <- rbind(
  phecodex_map  %>% select(phenoname, ukb_phenocode) ,
  ukb_aou_map %>% select(phenoname = aou_phenocode, ukb_phenocode)
) %>%
  merge(., v8_pheno_info %>% select(phenoname, description, category) %>% distinct(), by = 'phenoname', all.x = T)

table(aou_v8_mapped_ukb$category)
write_tsv(aou_v8_mapped_ukb, '~/Dropbox (Partners HealthCare)/aou_v8/data/aou_ukb_map/aou_ukb_matched_all_phenotype_v8.csv')
system("gsutil cp '/Users/wlu/Partners HealthCare Dropbox/Wenhan Lu/aou_v8/data/aou_ukb_map/aou_ukb_matched_all_phenotype_v8.csv' gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_ukb_matched_all_phenotype_v8.csv")
phenotype_maps <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/aou_ukb_map/aou_ukb_matched_all_phenotype_v8.csv')


aou_v8_mapped_pan_ukb <- rbind(
  phecode_map  %>% select(phenoname, ukb_phenocode = phecode) ,
  ukb_aou_map %>% select(phenoname = aou_phenocode, ukb_phenocode)
) %>%
  merge(., v8_pheno_info %>% select(phenoname, description, category) %>% distinct(), by = 'phenoname', all.x = T)
table(aou_v8_mapped_pan_ukb$category)
write_tsv(aou_v8_mapped_pan_ukb, '~/Dropbox (Partners HealthCare)/aou_v8/data/aou_ukb_map/aou_pan_ukb_matched_all_phenotype_v8.csv')
system("gsutil cp '/Users/wlu/Partners HealthCare Dropbox/Wenhan Lu/aou_v8/data/aou_ukb_map/aou_pan_ukb_matched_all_phenotype_v8.csv' gs://aou_wlu/v8_analysis/aou_ukb_meta/aou_pan_ukb_matched_all_phenotype_v8.csv")
phenotype_maps <- read_tsv('~/Dropbox (Partners HealthCare)/aou_v8/data/aou_ukb_map/aou_pan_ukb_matched_all_phenotype_v8.csv')
