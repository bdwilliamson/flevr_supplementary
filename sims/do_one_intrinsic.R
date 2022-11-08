## run the simulation once
#' @param indx the monte-carlo ID
#' @param n the sample size
#' @param p the number of features
#' @param family the type of outcome data to generate (e.g.,
#' continuous ["gaussian"] or binary ["binomial"])
#' @param link the link function
#' @param linear the type of linear predictor
#' (linear, e.g., X\beta; or nonlinear)
#' @param x_dist the type of x distribution
#' (i.e., "normal" for MVN or "nonnormal" for complex covariates)
#' @param outcome_regression whether or not the outcome regression should involve
#'   many functions ("many", the default) or a single, hard function ("hard")
#' @param rho the correlation between all features
#' @param rho2 the correlation between features in the active set
#' @param corr_type the type of correlation: "simple" (same rho used for all variables) or "complex" (negatively correlated active set)
#' @param beta the beta_0 parameter
#' @param missing_perc the percentage of missing data
#' (maximum column-wise, on average)
#' @param missing_type if missing_perc > 0, then MAR nested or MAR complex
#' @param estimator_type the type of estimator to run
#' (i.e., "lasso", "SL", "SPVIM", or "SAGE+")
#' @param extra_layer the type of added layer (e.g., "SS" or "KF")
#' @param learners the candidate learners (a list, with screens) to pass to SL
#' @param spvim_learners the candidate learners to pass to SPVIM
#' @param uni_learners the candidate learners to pass to SPVIM
#' for univariate regressions
#' @param b the number of bootstrap replicates
#' @param M the number of multiple imputations
#' @param alpha the per-comparison error rate
#' @param thresh the threshold for SL variable selection
#' (without SS or knockoffs)
#' @param rank the rank-based threshold for SL variable selection
#' (without SS or knockoffs)
#' @param ss_pfer the SS parameter controlling the per-family error rate
#' (expected number of false positives)
#' @param ss_pi the SS parameter controlling the threshold (specify either this or \code{ss_q})
#' @param ss_q the SS parameter controlling the number of variables considered (specify either this
#'   or \code{ss_pi})
#' @param kf_fdr the knockoff false discovery rate
#' @param pfp_q the parameter controlling the proportion of false positives
#' (for SPVIM)
#' @param gfwer_k the parameter controlling the number of type I
#' errors we can make (for SPVIM)
#' @param data_only only 'estimate' the true value on a large dataset
#' @param active_set the active set
#' @param impute_stability_threshold the threshold for the proportion of imputed datasets where a variable must be selected
do_one <- function(indx = 1, n = 100, p = 15, family = "binomial", link = NULL,
                   linear = TRUE, x_dist = "normal", outcome_regression = "many",
                   rho = 0, rho2 = 0, corr_type = "simple", beta = rep(0, p),
                   missing_perc = 0, missing_type = "nested",
                   estimator_type = "SL", extra_layer = "none",
                   learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                   spvim_learners = c("SL.ranger", "SL.glmnet", "SL.mean"),
                   uni_learners = "SL.polymars", b = 100, M = 1,
                   alpha = 0.05, thresh = 0.2, rank = 7,
                   ss_pfer = p * alpha, ss_pi = 0.75, ss_q = 7,
                   kf_fdr = 0.5,
                   pfp_q = 0.2, gfwer_k = 5, data_only = FALSE, active_set = 1:6,
                   impute_stability_threshold = 0.7) {
  # generate data
  if (linear) {
    data_gen <- gen_data_lm
  } else {
    data_gen <- gen_data_nlm
  }
  if (nchar(link) > 0) {
    family <- paste0(family, "-", link)
  }
  # get true performance
  test_dat <- data_gen(n = 5e4, p = p, rho = rho, rho2 = rho2, beta = matrix(beta), family = family,
                       x_dist = x_dist, y_dist = outcome_regression, corr_type = corr_type,
                       active_set = active_set)
  true_perf <- get_true_performance(test_dat = test_dat, linear = linear,
                                    family = family, beta = beta, x_dist = x_dist,
                                    y_dist = outcome_regression)
  if (data_only) {
    return(tibble::tibble(mc_id = indx, family = family, linear = linear,
                          x_dist = x_dist, missing_perc = missing_perc,
                          missing_type = missing_type, b = b, m_tot = M,
                          thresh = thresh, rank_thr = rank, n = n, p = p,
                          true_perf = true_perf))
  } 
  dat <- data_gen(n = n, p = p, rho = rho, rho2 = rho2, beta = matrix(beta), family = family,
                  x_dist = x_dist, y_dist = outcome_regression, corr_type = corr_type,
                  active_set = active_set)
  full_dat <- dat
  # if the percentage of missing data > 0, then make data missing
  if (missing_perc > 0) {
    missing_pattern <- make_missing_pattern(type = missing_type,
                                            p = p)
    # make data missing; may need to run this several times
    ampute_error <- TRUE
    while(ampute_error) {
      dat <- try(expr = mice::ampute(full_dat, prop = missing_perc,
                                     patterns = missing_pattern, mech = "MAR",
                          freq = mice::ampute.default.freq(missing_pattern),
                          weights = mice::ampute.default.weights(missing_pattern,
                                                                 "MAR"),
                          std = FALSE, bycases = TRUE) %>%
                   magrittr::use_series("amp") %>%
                   tibble::as_tibble(),
                 silent = TRUE)
      if (class(dat)[1] != "try-error") {
        ampute_error <- FALSE
      }
    }
    # multiply impute
    mice_preds <- mice::quickpred(data = dat, include = "y", mincor = 0.25)
    mice_impute <- mice::mice(data = dat, m = M, method = "pmm",
                              predictorMatrix = mice_preds,
                              where = is.na(dat), printFlag = FALSE)
    # create a list of lists: each element of imputed data is a list
    # with x and y
    imputed_data <- mice::complete(mice_impute, action = "long") %>%
      rename(imp = .imp, id = .id)
    # if using Long + Johnson, bootstrap first
    if (grepl("LJ", estimator_type) | grepl("BI-BL", estimator_type)) {
      boot_data_list <- generate_bootstrapped_data(data = dat, B = b)
      boot_imputed_data <- lapply(boot_data_list, function(dat) {
        mice_preds <- mice::quickpred(data = dat, include = "y", mincor = 0.25)
        mice_impute <- mice::mice(data = dat, m = 1, method = "pmm",
                                  predictorMatrix = mice_preds,
                                  where = is.na(dat), printFlag = FALSE)
        # create a list of lists: each element of imputed data is a list
        # with x and y
        imputed_dat <- mice::complete(mice_impute, action = "long") %>%
          rename(imp = .imp, id = .id)
        imputed_dat
      })
    } 
  } else { # all deltas are zero
    # no need to multiply impute
    full_dat <- dat
    imputed_data <- dat %>%
      mutate(id = 1:nrow(dat), imp = 1)
  }
  # set up return tibble
  selected_sets <- NULL
  est_vim <- NULL
  fam_for_ests <- switch((grepl("binomial", family)) + 1,
                         "gaussian", "binomial")
  # get a set of selected variables
  opts <- get_sl_opts(fam_for_ests)
  select_data <- switch(as.numeric(grepl("LJ", estimator_type) | grepl("BI-BL", estimator_type)) + 1, imputed_data, boot_imputed_data)
  selection_results <- select_from_all(
    imputed_training_data = select_data, outcome_name = "y", M = M, 
    variable_selection_procedure = paste0(estimator_type, "_", extra_layer),
    impute_stability_threshold = impute_stability_threshold, alpha = alpha,
    sl_rank_threshold = rank, ss_pfer = ss_pfer, ss_pi = ss_pi, 
    kf_fdr = kf_fdr, pfp_q = pfp_q, gfwer_k = gfwer_k, ss_b = b,
    SL.library = learners, spvim_library = spvim_learners, 
    univariate_SL.library = uni_learners, method = "method.CC_nloglik2", 
    family = opts$fam, method = opts$method, 
    cvControl = list(V = 5, stratifyCV = TRUE)
  )
  selected_set <- selection_results$pooled
  # get performance of selected set
  if (!is.list(selected_set)) {
    est_perf_tib <- data.table::data.table(imp = 1:M, type = estimator_type, 
                                           perf = NA_real_, fit_out = NA_character_, 
                                           sens = NA_real_, spec = NA_real_, fdr = NA_real_)
    selected_var_str <- paste(which(selected_set == 1), collapse = "; ")
  } else {
    est_perf_tib <- data.table::as.data.table(expand.grid(imp = 1:M, type = c("gFWER", "PFP", "FDR"))) %>%
      mutate(perf = NA_real_, fit_out = NA_character_, sens = NA_real_, spec = NA_real_, fdr = NA_real_)
    selected_var_str <- lapply(
      selected_set, function(z) paste(which(z == 1), collapse = "; ")
    )
    names(selected_set) <- c("gFWER", "PFP", "FDR")
  }
  for (m in seq_len(M)) {
    x <- imputed_data %>% filter(imp == m) %>% select(-y, -imp, -id)
    y <- imputed_data %>% filter(imp == m) %>% pull(y)
    test_x <- test_dat %>% select(-y)
    test_y <- test_dat %>% pull(y)
    no_screen_lib <- unlist(lapply(learners,
                                   function(x) x[!grepl("screen", x) &
                                                   !grepl("All", x)]))
    if (!is.list(selected_set)) {
      selected_mod_perf <- get_selected_set_perf(
        selected_set = selected_set, x = x, y = y,
        test_x = test_x, test_y = test_y,
        estimator_type = estimator_type,
        learner_lib = no_screen_lib, family = fam_for_ests
      )
      est_perf_tib[m, ]$perf <- selected_mod_perf$perf
      est_perf_tib[m, ]$fit_out <- selected_mod_perf$fit_out
      est_perf_tib[m, ]$sens <- get_sensitivity(selected_var_str, active_set)
      est_perf_tib[m, ]$spec <- get_specificity(selected_var_str, active_set, p)
      est_perf_tib[m, ]$fdr <- sum((selected_set * c(rep(0, length(active_set)), rep(1, p - length(active_set)))) / sum(selected_set))
    } else {
      selected_mod_perf <- sapply(
        1:length(selected_set),
        function(s) {
          get_selected_set_perf(
            selected_set = selected_set[[s]], x = x, y = y, test_x = test_x,
            test_y = test_y, estimator_type = estimator_type,
            learner_lib = no_screen_lib, family = fam_for_ests
          )
        }, simplify = FALSE)
      est_perf_tib[imp == m, 'perf'] <- unlist(lapply(selected_mod_perf, function(l) l$perf))
      est_perf_tib[imp == m, 'fit_out'] <- unlist(lapply(selected_mod_perf, function(l) l$fit_out))
      est_perf_tib[imp == m, 'sens'] <- unlist(lapply(
        selected_var_str, get_sensitivity, true_active_set = active_set
      ))
      est_perf_tib[imp == m, 'spec'] <- unlist(lapply(
        selected_var_str, get_specificity, true_active_set = active_set, p = p
      ))
      est_perf_tib[imp == m, 'fdr'] <- unlist(lapply(
        selected_set, function(l) sum((l * c(rep(0, length(active_set)), rep(1, p - length(active_set)))) / sum(l))
      ))
    }
  }
  # combine performance using Rubin's rules
  est_perf <- est_perf_tib %>% 
    group_by(type) %>% 
    summarize(perf = mean(perf), fit_out = first(fit_out), 
              sens = mean(sens), spec = mean(spec), fdr = mean(fdr), .groups = "drop")
  miss_percs <- colMeans(is.na(dat))
  miss_summ <- summary(miss_percs)
  miss_summ_text <- sprintf(
    "Min. = %1$.3f; 1st Qu. = %2$.3f; Median = %3$.3f; Mean = %4$.3f; 3rd Qu. = %5$.3f; Max. = %6$.3f.",
    miss_summ[1], miss_summ[2], miss_summ[3],
    miss_summ[4], miss_summ[5], miss_summ[6]
  )
  this_est_type <- switch((is.list(selected_set)) + 1,
                          estimator_type,
                          paste(estimator_type, names(selected_set), sep = "-"))
  set_tib <- tibble::tibble(estimator_type = this_est_type,
                            extra_layer = extra_layer, m = m,
                            selected_vars = selected_var_str,
                            perf = est_perf$perf,
                            true_perf = true_perf,
                            sensitivity = est_perf$sens, specificity = est_perf$spec,
                            fdr = est_perf$fdr,
                            overall_miss = mean(miss_percs),
                            miss_summ = miss_summ_text,
                            final_fit_out = est_perf$fit_out)
  set_output <- set_tib %>% 
    tibble::add_column(mc_id = indx, family = family, linear = linear,
                       x_dist = x_dist, missing_perc = missing_perc,
                       missing_type = missing_type, b = b, m_tot = M,
                       thresh = thresh, rank_thr = rank, n = n, p = p,
                       .before = "estimator_type")
  vim_output <- selection_results$vim %>% 
    mutate(estimator_type = estimator_type, .before = "feature") %>% 
    tibble::add_column(mc_id = indx, family = family, linear = linear,
                       x_dist = x_dist, missing_perc = missing_perc,
                       missing_type = missing_type, b = b, m_tot = M,
                       thresh = thresh, rank_thr = rank, n = n, p = p,
                       .before = "estimator_type")
  list(selected = set_output, vim = vim_output)
}
