# compile results into a single file

# set up -----------------------------------------------------------------------
library("here")
library("argparse")
library("data.table")
library("tibble")
library("dplyr")

# get command-line args
parser <- ArgumentParser()
parser$add_argument("--M", type = "double", default = 10, 
                    help = "number of imputations")
parser$add_argument("--est-type",
                    default = "glm",
                    help = "estimation procedure to use")
parser$add_argument("--selection-type",
                    default = "lasso",
                    help = "variable selection procedure to use")
parser$add_argument("--algorithm", default = "1",
                    help = "which algorithm for pooling to use?")
parser$add_argument("--output-dir",
                    default = "results/data_analysis/",
                    help = "directory to send output to")
args <- parser$parse_args()

# read in all of the results for the given selection type + estimator + algorithm combo -----
output_dir <- paste0(args$output_dir, "output_", args$selection_type,
                     "_", args$est_type, "_", args$algorithm, "/",
                     "impute_", args$M, "/")
perf_files <- paste0(output_dir, list.files(output_dir, pattern = "perf"))
vim_files <- paste0(output_dir, list.files(output_dir, pattern = "vim"))

perf <- as_tibble(rbindlist(
    lapply(as.list(seq_len(length(perf_files))), function(l) {
        this_file <- perf_files[l]
        this_outcome <- ifelse(grepl("mucinous", this_file), "mucinous", "high_malignancy")
        covariates <- as.numeric(grepl(paste0(this_outcome, "_1"), this_file))
        biomarkers <- as.numeric(grepl("1_id", this_file))
        readRDS(this_file) %>% 
            mutate(outcome = this_outcome,
                   covariates = covariates, biomarkers = biomarkers,
                   .before = "selection_procedure")
    })
))
vim <- as_tibble(rbindlist(
    lapply(as.list(seq_len(length(vim_files))), function(l) {
        this_file <- vim_files[l]
        this_outcome <- ifelse(grepl("mucinous", this_file), "mucinous", "high_malignancy")
        covariates <- as.numeric(grepl(paste0(this_outcome, "_1"), this_file))
        biomarkers <- as.numeric(grepl("1_id", this_file))
        readRDS(this_file) %>% 
            mutate(outcome = this_outcome,
                   covariates = covariates, biomarkers = biomarkers,
                   .before = "selection_procedure")
    })
))

# save
common_name <- paste0(args$selection_type, "_", args$est_type, "_", args$algorithm, ".rds")
saveRDS(perf, file = here("results", "data_analysis", paste0("perf_", common_name)))
saveRDS(vim, file = here("results", "data_analysis", paste0("vim_", common_name)))
