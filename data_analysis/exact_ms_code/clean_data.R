# clean the initial dataset into an analysis dataset

# load required libraries ------------------------------------------------------
library("here")
library("tidyverse")
library("janitor")

source(here("code", "00_utils.R"))

# args
prop_miss <- 10
pairwise_corr <- 0.9

# read in the data -------------------------------------------------------------
raw_outcome_biomarker_data <- readr::read_csv(here("data", "data_merge_new.csv")) %>% 
  janitor::clean_names()
raw_covariate_data <- readxl::read_xlsx(here("data", "PCBA_combined_clinical_data_2020-03-26.xlsx"), 
                                        sheet = 1) %>% 
  janitor::clean_names()

# clean up the outcome/biomarker data ------------------------------------------
# interested in site (institution) and all biomarkers
# from the data for analysis 1 and/or analysis 2; these adjudicate the outcomes
# ignore warnings, these are due to creating numerics from the character vectors
outcome_biomarker_data <- raw_outcome_biomarker_data %>% 
  select(id, institution, analysis_1, analysis_2,
         starts_with("lab1_"), starts_with("lab2_"), starts_with("lab3_"),
         starts_with("lab4_"), starts_with("lab5_"), starts_with("lab6_"),
         starts_with("cea")) %>% 
  rename(mucinous = analysis_1, high_malignancy = analysis_2) %>% 
  mutate(lab4_glucose_score = case_when(
    lab4_glucose_score == "<20" ~ 10,
    lab4_glucose_score == "Not measurable" ~ NA_real_,
    !is.na(as.numeric(lab4_glucose_score)) ~ as.numeric(lab4_glucose_score)
  ), lab4_areg_score = case_when(
    lab4_areg_score == ">4000" ~ 8000,
    lab4_areg_score == "Not measurable" ~ NA_real_,
    lab4_areg_score == "Not Detectable" ~ 0.5,
    !is.na(as.numeric(lab4_areg_score)) ~ as.numeric(lab4_areg_score)
  ))


# merge datasets together ------------------------------------------------------
# remove calls that I can't verify
merged_clean_data <- outcome_biomarker_data %>% 
  select(-lab3_muc3ac_mucinous_call, -lab3_muc5ac_mucinous_call,
         -lab3_mucinous_both_call) %>% 
  dplyr::left_join(screened_covariate_data, by = c("id", "institution"))

# save off the clean data for both outcomes
saveRDS(merged_clean_data %>% 
          select(-high_malignancy, -id), 
        here("analysis_data", "objective_1_data.rds"))
saveRDS(merged_clean_data %>% 
          select(-mucinous, -id), 
        here("analysis_data", "objective_2_data.rds"))
