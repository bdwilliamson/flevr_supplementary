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
parser$add_argument("--nreps-total", type = "double", default = 100,
                    help = "total number of replicates")
parser$add_argument("--nreps-per-job", type = "double", default = 5,
                    help = "total number of replicates per job")
parser$add_argument("--K", type = "double",
                    default = 5,
                    help = "number of cross-fitting folds")
parser$add_argument("--M", type = "double", default = 10,
                    help = "number of imputations")
parser$add_argument("--est-type",
                    default = "SL",
                    help = "estimation procedure to use")
parser$add_argument("--selection-type",
                    default = "SPVIM-RR",
                    help = "variable selection procedure to use")
parser$add_argument("--algorithm", default = "1",
                    help = "which algorithm for pooling to use?")
parser$add_argument("--output-dir",
                    default = "results/",
                    help = "directory to send output to")
args <- parser$parse_args()

# set up args, get functions ---------------------------------------------------
# pull in job id; edit the following line if using a different system than Slurm
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(job_id)) {
  job_id <- 1
  code_dir <- switch(
    grepl("code", here()) + 1, "code", ""
  )
  data_folder <- switch(
    grepl("code", here()) + 1, "analysis_data", "../analysis_data"
  )
  args$output_dir <- here("results", paste0("output_",
                                            args$selection_type, "_",
                                            args$est_type, "_",
                                            args$algorithm))
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
nreps_per_combo <- args$nreps_total / args$nreps_per_job
covariates <- c(0)
biomarkers <- c(1)
all_outcomes <- c("mucinous", "high_malignancy")
all_analyses <- tibble::as_tibble(
  expand.grid(mc_id = 1:nreps_per_combo,
              outcome = all_outcomes, covariates = covariates,
              biomarkers = biomarkers, stringsAsFactors = FALSE)
)
all_analyses <- all_analyses %>%
  filter(!(covariates == 0 & biomarkers == 0))
this_analysis <- all_analyses[job_id, ]

biomarker_covariate_txt <- switch(this_analysis$biomarkers + 1, "risk factors only,",
                                  paste0("biomarkers",
                                         switch(this_analysis$covariates + 1,
                                                " only,", " and other risk factors,")))
cat("Running analysis procedure", args$selection_type, "+",
    args$est_type, "on outcome", paste0("'", this_analysis$outcome, "'"),
    "using", biomarker_covariate_txt,
    "with", args$nreps_per_job, "monte-carlo replications and",
    paste0(args$K, "-fold cross-fitting.\n"))

output_dir <- paste0(args$output_dir, "/impute_", args$M, "/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
common_filename <- paste0(this_analysis$outcome, "_", this_analysis$covariates,
                          "_", this_analysis$biomarkers,
                          "_id_", this_analysis$mc_id, ".rds")
covariate_filename <- paste0(this_analysis$outcome, "_1_0_id_", this_analysis$mc_id, ".rds")

# read in the data -------------------------------------------------------------
dataset_name <- switch(
  grepl("mucinous", this_analysis$outcome) + 1, "objective_2_data.rds",
  "objective_1_data.rds"
)
dataset <- readRDS(
  here(data_folder, dataset_name)
)
# subset to only biomarkers if specified by the analysis
if (this_analysis$covariates == 0) {
  dataset <- dataset %>%
    select(institution, !!this_analysis$outcome, starts_with("jhu"),
           starts_with("ucsf"), starts_with("vandel"), starts_with("stanford"),
           starts_with("upmc"), starts_with("washu"), starts_with("cea"))
}
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

if (grepl("screen", args$est_type)) {
  screen_lib <- c("screen.glmnet", "screen.corRank.p", "All")
  learner_lib <- lapply(as.list(learner_lib),
                        function(x) c(x, screen_lib))
  final_learner_lib <- lapply(as.list(final_learner_lib),
                              function(x) c(x, screen_lib))
}

cat("Estimating performance using Monte-Carlo replicates and cross-fitting\n")

# set up static args
impute_threshold <- 0.7
target_alpha <- 0.04
p <- ncol(dataset) - 1
ss_pi <- 0.9
ss_pfer <- p * target_alpha
sl_rank_thresh <- ceiling(sqrt(ss_pfer * (2 * ss_pi - 1) * p)) + 5
gfwer_k <- 5
pfp_q <- 0.8
kf_fdr <- 0.2
ss_b <- 100

# compute the random number seed
all_algos <- c("SPVIM-RR", "lasso-LJ", "lasso-BI-BL")
reg_procs <- c("SL", "glm", "glm")
all_procedures <- tibble::tibble(
  selection_procedure = rep(all_algos, each = 2),
  estimation_procedure = rep(reg_procs, each = 2),
  cv_alg = rep(c(1, 2), length(all_algos))
) %>% mutate(procedure = row_number())
# note that seed runs from 1:(nreps_per_combo * 14 [= number unique procedures])
# want a consistent set of seeds for each unique procedure
which_procedure <- all_procedures %>%
  filter(selection_procedure == args$selection_type, estimation_procedure == args$est_type,
         cv_alg == as.numeric(args$algorithm)) %>% pull(procedure)
print(which_procedure)
print(job_id)
seed <- which_procedure * 1e3 + job_id
cat("Seed:", seed, "\n")
set.seed(seed)
# actually run the analysis
system.time(
  output <- lapply(
    1:args$nreps_per_job, function(i) {
      run_analysis_once(
        mc_id = i + args$nreps_per_job * (this_analysis$mc_id - 1),
        dataset = dataset, preimputed_data = preimputed_data,
        outcome_name = this_analysis$outcome,
        variable_selection_procedure = args$selection_type,
        prediction_procedure = args$est_type,
        covariates = (this_analysis$covariates == 1), biomarkers = (this_analysis$biomarkers == 1),
        K = args$K, M = args$M, impute_stability_threshold = impute_threshold,
        alpha = target_alpha, p = p, sl_rank_threshold = sl_rank_thresh, ss_b = ss_b,
        ss_pfer = p * target_alpha, ss_pi = ss_pi, ss_q = sl_rank_thresh,
        kf_fdr = kf_fdr, pfp_q = pfp_q, gfwer_k = sl_rank_thresh,
        algorithm = args$algorithm, SL.library = final_learner_lib,
        spvim_library = spvim_lib, univariate_SL.library = uni_learner_lib,
        method = "method.CC_nloglik2", family = binomial(),
        cvControl = list(V = 5, stratifyCV = TRUE), missing_threshold = 0.1
      )
    }
  )
)
# break off performance, estimated VIM, imputed data
perf <- as_tibble(rbindlist(lapply(output, function(l) l$perf)))
vim <- as_tibble(rbindlist(lapply(output, function(l) l$vim)))
imputed_data <- as_tibble(rbindlist(lapply(output, function(l) {
  as_tibble(rbindlist(l$data))
})))
# save it
saveRDS(perf, paste0(output_dir, "perf_", common_filename))
saveRDS(vim, paste0(output_dir, "vim_", common_filename))
saveRDS(imputed_data, paste0(output_dir, "imputed_data_", common_filename))
print("Analysis complete!")
