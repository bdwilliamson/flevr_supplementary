# run a single round of the analysis
# @param mc_id the monte-carlo ID
# @param dataset the dataset to use
# @param preimputed_data pre-imputed covariate data (only used if covariates and biomarkers are being analyzed)
# @param outcome_name the name of the outcome
# @param variable_selection_procedure the variable selection procedure (e.g., lasso)
# @param prediction_procedure the prediction procedure (e.g., glm)
# @param covariates logical, are we using covariates?
# @param biomarkers logical, are we using biomarkers?
# @param K the number of cross-fitting folds
# @param M the number of imputations
# @param impute_stability_threshold the proportion of imputed datasets in which
#   a variable must be selected to be part of the final set of variables
# @param alpha the desired per-comparison error rate
# @param sl_rank_threshold the rank-based threshold for SL variable selection;
#   by default,
# @param ss_pfer the SS parameter controlling the per-family error rate (expected num. false positives)
# @param ss_pi the SS parameter controlling the threshold (specify either this or \code{ss_q})
# @param ss_q the SS parameter controlling the number of variables considered (specify either this or \code{ss_pi})
# @param kf_fdr the desired knockoff FDR
# @param pfp_q the desired proportion of false positives (for intrinsic selection)
# @param gfwer_k the desired number of generalized family-wise errors we can make (for intrinsic selection)
# @param algorithm the algorithm to use; in Alg 1, we pool before training prediction algos;
#   in Alg 2, we don't.
# @param missing_threshold the proportion of missing data to allow CV-AUC computation
# @param ss_b the number of bootstrap replicates for stability selection-type algorithms
# @param ... other arguments to pass to the variable selection and prediction procedures
run_analysis_once <- function(mc_id = 1, dataset, preimputed_data = NULL,
                              outcome_name = "mucinous",
                              variable_selection_procedure = "lasso",
                              prediction_procedure = "glm",
                              covariates = TRUE, biomarkers = TRUE,
                              K = 5, M = 10, impute_stability_threshold = 0.7,
                              alpha = 0.04, p = ncol(dataset) - 1,
                              ss_b = 100, ss_pfer = p * alpha, ss_pi = 0.9,
                              sl_rank_threshold = ceiling(sqrt(ss_pfer * (2 * ss_pi - 1) * p)) + 1,
                              ss_q = sl_rank_threshold, kf_fdr = 0.2,
                              pfp_q = 0.8, gfwer_k = ss_q,
                              algorithm = "1", missing_threshold = 0.1, 
                              ...) {
  print(mc_id)
  arg_lst <- list(...)
  # set up the data
  y <- dataset %>% pull(!!outcome_name)
  # create folds for cross-fitting
  if (sum(is.na(y)) > 0) {
    cross_fitting_folds <- vector("numeric", length(y))
    cross_fitting_folds[!is.na(y)] <- vimp::make_folds(y = y[!is.na(y)], V = K, stratified = TRUE)
    cross_fitting_folds[is.na(y)] <- sample(1:K, size = sum(is.na(y)))
  } else {
    cross_fitting_folds <- vimp::make_folds(y = y, V = K, stratified = TRUE)
  }
  bootstrap_first <- grepl("LJ", variable_selection_procedure) | 
    grepl("BI-BL", variable_selection_procedure)
  # set up lists to collect results for returning
  perf_lst <- vector("list", length = K)
  vim_lst <- vector("list", length = K)
  imputed_data_lst <- vector("list", length = K)
  impute_vars <- names(dataset)[!grepl("institution", names(dataset))]
  dataset <- dataset %>% mutate(id = 1:nrow(dataset), .before = "institution")
  # do cross-fitting
  for (k in seq_len(K)) {
    print(k)
    # get imputed training and testing data
    all_imputed_data <- impute_train_and_test(dataset = dataset, preimputed_data = preimputed_data,
                                              mc_id = mc_id, cross_fitting_folds = cross_fitting_folds,
                                              k = k, M = M, covariates = covariates,
                                              biomarkers = biomarkers,
                                              no_impute_vars = switch(as.numeric(covariates & biomarkers) + 1,
                                                                      c("id", "institution"),
                                                                      c("id", "mc_id", "k", "institution")),
                                              use_cea = !grepl("high_malignancy", outcome_name),
                                              bootstrap_first = bootstrap_first,
                                              B = ss_b)
    imputed_training_data <- all_imputed_data$train
    imputed_testing_data <- all_imputed_data$test
    imputed_data_lst[[k]] <- all_imputed_data$all
    boot_imputed_training_data <- all_imputed_data$train_boot
    # do variable selection:
    #   if using SPVIM with Rubin's rules, estimate VIMs, pool, and do one variable selection
    #   otherwise, do variable selection based on each imputed dataset separately
    train_for_selecting <- switch(as.numeric(bootstrap_first) + 1, imputed_training_data, 
                                  boot_imputed_training_data)
    selection_results <- select_from_all(imputed_training_data = train_for_selecting,
                                         outcome_name = outcome_name, M = M,
                                         variable_selection_procedure = variable_selection_procedure,
                                         impute_stability_threshold = impute_stability_threshold,
                                         alpha = alpha, sl_rank_threshold = sl_rank_threshold,
                                         ss_pfer = ss_pfer, ss_pi = ss_pi, ss_q = ss_q,
                                         kf_fdr = kf_fdr, pfp_q = pfp_q, gfwer_k = gfwer_k,
                                         ss_b = ss_b, ...)
    pooled_set <- selection_results$pooled
    selected_vars_lst <- selection_results$all
    vim_lst[[k]] <- selection_results$vim
    # estimate prediction performance:
    #   if Alg 1, use the pooled set
    #   if Alg 2, use each individual set (not available for SPVIM-RR)
    if (algorithm == "1") {
        variable_sets <- rep(list(pooled_set), M)
    } else {
        variable_sets <- selected_vars_lst
    }
    performance_lst <- prediction_performance(M = M, train = imputed_training_data,
                                              test = imputed_testing_data,
                                              variable_sets = variable_sets,
                                              n = length(y),
                                              outcome_name = outcome_name,
                                              prediction_procedure = prediction_procedure,
                                              arg_lst = arg_lst,
                                              missing_threshold = missing_threshold)
    # combine performance using Rubin's rules
    perf_lst[[k]] <- pool_all_ests(performance_lst = performance_lst, variable_selection_procedure = variable_selection_procedure, k = k, n = nrow(dataset))
  }
  # collapse
  if (grepl("SPVIM", variable_selection_procedure)) {
    grouped_perf <- as_tibble(rbindlist(perf_lst)) %>%
      group_by(n, type)
    procedure <- paste0(variable_selection_procedure, " + ", c("FDR", "gFWER", "PFP"))
  } else {
    grouped_perf <- as_tibble(rbindlist(perf_lst)) %>%
      group_by(n)
    procedure <- variable_selection_procedure
  }
  all_perf <- grouped_perf %>%
    summarize(est = mean(est, na.rm = TRUE), var = mean(var, na.rm = TRUE), .groups = "drop") %>%
    bind_cols(tibble::tibble(mc_id = mc_id, selection_procedure = procedure,
              pred_alg = prediction_procedure, algo = algorithm)) %>%
    select(mc_id, selection_procedure, pred_alg, algo, n, est, var)
  all_vim <- as_tibble(rbindlist(vim_lst)) %>%
    bind_cols(tibble::tibble(mc_id = mc_id, selection_procedure = variable_selection_procedure)) %>%
    select(mc_id, selection_procedure, everything())
  list(perf = all_perf, vim = all_vim, data = imputed_data_lst)
}
