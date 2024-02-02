#!/usr/local/bin/Rscript

# run the analysis for a given estimator and dataset

# load required packages -------------------------------------------------------
library("here") # navigation
library("argparse") # command-line args
library("tidyverse") # nice functions
library("data.table") # nice functions
library("mice") # imputation
library("SuperLearner") # algos
library("glmnet")
library("ranger")
library("xgboost")
library("kernlab")
lib_paths <- .libPaths()
# variable importance
library("vimp")
# variable selection
library("flevr")
library("stabs")
library("knockoff")
library("parallel")
# cross-validated AUC
library("cvAUC")

# get command-line args
parser <- ArgumentParser()
parser$add_argument("--M", type = "double", default = 10,
                    help = "number of imputations")
parser$add_argument("--selection-type",
                    default = "lasso",
                    help = "variable selection procedure to use")
parser$add_argument("--output-dir",
                    default = "results",
                    help = "directory to send output to")
args <- parser$parse_args()

# set up args, get functions ---------------------------------------------------
# pull in job id; edit the following line if not using Slurm
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(job_id)) {
    job_id <- 1
    code_dir <- switch(
        grepl("code", here()) + 1, "code", ""
    )
    data_folder <- switch(
        grepl("code", here()) + 1, "analysis_data", "../analysis_data"
    )
    args$output_dir <- here("results", "data_analysis", 
                            paste0("output_", args$selection_type))
} else {
    code_dir <- ""
    data_folder <- "analysis_data"
}

# read in helper functions
source(here(code_dir, "01_impute_missing_data.R"))
source(here(code_dir, "02_variable_selection.R"))
source(here(code_dir, "03_prediction_performance.R"))
source(here(code_dir, "04_run_analysis.R"))
source(here(code_dir, "00_utils.R"))

# set up all possible analyses
covariates <- c(0)
biomarkers <- c(1)
all_outcomes <- c("mucinous", "high_malignancy")
all_analyses <- tibble::as_tibble(
    expand.grid(outcome = all_outcomes, covariates = covariates,
                biomarkers = biomarkers, stringsAsFactors = FALSE)
)
all_analyses <- all_analyses %>%
    filter(!(covariates == 0 & biomarkers == 0))
this_analysis <- all_analyses[job_id, ]

biomarker_covariate_txt <- switch(this_analysis$biomarkers + 1, "risk factors only.\n",
                                  paste0("biomarkers",
                                         switch(this_analysis$covariates + 1,
                                                " only.\n", " and other risk factors.\n")))
cat("Running analysis procedure", args$selection_type, "on outcome", paste0("'", this_analysis$outcome, "'"),
    "using", biomarker_covariate_txt)

output_dir <- paste0(args$output_dir, "/impute_", args$M, "/")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# read in the data -------------------------------------------------------------
dataset_name <- switch(
    grepl("mucinous", this_analysis$outcome) + 1, "objective_2_data.rds",
    "objective_1_data.rds"
)
dataset_init <- readRDS(
    here(data_folder, dataset_name)
)
dataset <- dataset_init %>%
    mutate(id = 1:nrow(dataset_init), .before = "institution")
analysis_dataset <- dataset
# subset to only biomarkers if specified by the analysis
if (this_analysis$covariates == 0) {
    analysis_dataset <- dataset %>%
        select(id, institution, !!this_analysis$outcome, starts_with("jhu"),
               starts_with("ucsf"), starts_with("vandel"), starts_with("stanford"),
               starts_with("upmc"), starts_with("washu"), starts_with("cea"))
}

common_filename <- paste0(this_analysis$outcome, "_", this_analysis$covariates,
                          "_", this_analysis$biomarkers,
                          "_final.rds")

covariate_filename <- paste0(this_analysis$outcome, "_1_0_final.rds")
# read in the imputed data, if covariates and biomarkers are being used
preimputed_data <- NULL

# run the analysis -------------------------------------------------------------
# set up the SL library
xgb_tune_params <- list(max_depth = 4,
                        shrinkage = 0.1,
                        nrounds = c(100, 500, 2000))
xgb_learners <- create.Learner("SL.xgboost_new", tune = xgb_tune_params,
                               detailed_names = TRUE, name_prefix = "xgb")
ranger_tune_params = list(
    max.depth = c(1, 10, 20, 30, 100)
)
rf_learners <- create.Learner("SL.ranger.imp", tune = ranger_tune_params,
                              detailed_names = TRUE, name_prefix = "rf")
learner_lib <- c(paste0("SL.glmnet.", c(0, 25, 50, 75)),
                 xgb_learners$names,
                 rf_learners$names, "SL.ranger.imp")
final_learner_lib <- learner_lib
spvim_lib <- lapply(as.list(xgb_learners$names), function(x) c(x, "screen.corRank.p"))
uni_learner_lib <- "my_SL.polymars"

# set up static args
impute_threshold <- 0.7
target_alpha <- 0.04
p <- ncol(analysis_dataset) - 1
ss_pi <- 0.9
ss_pfer <- p * target_alpha
sl_rank_thresh <- ceiling(sqrt(ss_pfer * (2 * ss_pi - 1) * p)) + 5
gfwer_k <- 5
pfp_q <- 0.8
kf_fdr <- 0.2
ss_b <- 100

# compute the random number seed
all_algos <- c("lasso-LJ", "lasso-BI-BL", "SPVIM-RR")
# note that seed runs from 1:(nreps_per_combo * 14 [= number unique procedures])
# want a consistent set of seeds for each unique procedure
seed <- 1234
cat("Seed:", seed, "\n")
set.seed(seed)

cat("Imputing data\n")
# actually run the analysis
bootstrap_first <- grepl("LJ", args$selection_type) | grepl("BI-BL", args$selection_type)
system.time(
    if (this_analysis$covariates & this_analysis$biomarkers) {
        covariate_nms <- names(preimputed_data %>% select(-imp, -id))
        biomarkers_and_imputed_covars <- lapply(as.list(unique(preimputed_data$imp)), function(i) {
            left_join(analysis_dataset %>% select(-all_of(covariate_nms)),
                      preimputed_data %>% filter(imp == !!i), by = "id")
        })
        imputed_datasets <- as_tibble(rbindlist(lapply(as.list(seq_len(args$M)), function(m) {
            these_data <- biomarkers_and_imputed_covars[[m]] %>%
                filter(imp == !!m)
            imputed_init <- impute_missing_data(dataset = these_data %>% select(-imp), M = 1,
                                                covariates_only = FALSE,
                                                no_impute_vars = c("id", "institution"),
                                                use_cea = !grepl("high_malignancy", this_analysis$outcome))
            imputed_init %>% select(-imp) %>% mutate(imp = these_data$imp)
        })))
    } else {
        if (bootstrap_first) {
          boot_data_list <- generate_bootstrapped_data(data = analysis_dataset, B = ss_b)
          boot_imputed_data <- lapply(boot_data_list, function(dat) {
            impute_missing_data(dataset = dat, M = 1, 
                                covariates_only = this_analysis$covariates & !this_analysis$biomarkers, 
                                no_impute_vars = c("id", "institution"),
                                use_cea = !grepl("high_malignancy", this_analysis$outcome))
          })
        } else {
          boot_imputed_data <- NULL
        }
        imputed_datasets <- impute_missing_data(dataset = analysis_dataset, M = args$M,
                                                covariates_only = (this_analysis$covariates == 1 & this_analysis$biomarkers == 0),
                                                no_impute_vars = c("id", "institution"),
                                                use_cea = !grepl("high_malignancy", this_analysis$outcome))
    }
)
saveRDS(imputed_datasets, paste0(output_dir, "imputed_data", common_filename))
saveRDS(boot_imputed_data, paste0(output_dir, "boot_imputed_data", common_filename))

dat_for_selecting <- switch(as.numeric(bootstrap_first) + 1, imputed_datasets, 
                            boot_imputed_data)
cat("Selecting variables\n")
system.time(
    selection_results <- select_from_all(imputed_training_data = dat_for_selecting,
                                         outcome_name = this_analysis$outcome, M = args$M,
                                         variable_selection_procedure = args$selection_type,
                                         impute_stability_threshold = impute_threshold,
                                         alpha = target_alpha, p = p, sl_rank_threshold = sl_rank_thresh, ss_b = ss_b,
                                         ss_pfer = p * target_alpha, ss_pi = ss_pi, ss_q = sl_rank_thresh,
                                         kf_fdr = kf_fdr, pfp_q = pfp_q, gfwer_k = sl_rank_thresh,
                                         algorithm = args$algorithm, SL.library = final_learner_lib,
                                         spvim_library = spvim_lib, univariate_SL.library = uni_learner_lib,
                                         method = "method.CC_nloglik2", family = binomial(),
                                         cvControl = list(V = 5, stratifyCV = TRUE))
)
pooled_set <- selection_results$pooled
vim <- selection_results$vim
saveRDS(vim, paste0(output_dir, "vim_", common_filename))

# final set of variables and the estimated VIMs --------------------------------
if (is.list(pooled_set)) {
    final_selected_vars_lst <- lapply(pooled_set, function(set) {
        names(analysis_dataset %>% select(-id, -!!this_analysis$outcome))[set]
    })
    procedures <- paste0(args$selection_type, " + ", c("gFWER", "PFP", "FDR"))
    final_selected_vars <- as_tibble(rbindlist(lapply(as.list(1:3), function(l) {
        tibble(procedure = procedures[l],
               selected = final_selected_vars_lst[[l]])
    })))
} else {
    final_selected_vars_lst <- names(analysis_dataset %>%
                                     select(-id, -!!this_analysis$outcome))[pooled_set == 1]
    final_selected_vars <- tibble(procedure = args$selection_type,
                                  selected = final_selected_vars_lst)
}
saveRDS(final_selected_vars, paste0(output_dir, "selected-vars_", common_filename))
print("Analysis complete!")
