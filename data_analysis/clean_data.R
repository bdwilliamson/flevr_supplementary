# clean the initial dataset into an analysis dataset

# load required libraries ------------------------------------------------------
library("flevr")

# read in the data -------------------------------------------------------------
data(biomarkers)

# save off the clean data for both outcomes
saveRDS(biomarkers %>% 
          select(-high_malignancy, -id), 
        here("analysis_data", "objective_1_data.rds"))
saveRDS(biomarkers %>% 
          select(-mucinous, -id), 
        here("analysis_data", "objective_2_data.rds"))
