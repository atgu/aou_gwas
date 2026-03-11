source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
if (!require("ggpattern")) install.packages("ggpattern")


pop_pattern = c('afr' = 'white',
                'amr' = 'grey90',
                'eas' = 'white',
                'eur' = 'grey90',
                'mid' = 'white',
                'sas'= 'grey90',
                'aou_meta'= 'white',
                'ukb_meta' = 'grey90',
                'all_meta'= 'white')

full_new_hits <- obtain_novelty_data_for_plot()

plot_novelty(data= full_new_hits, category_name='lab_measurement', pop_pattern=pop_pattern, indexes = 1:9)
plot_novelty(data= full_new_hits, category_name='physical_measurement', pop_pattern=pop_pattern, indexes = 1:9)
