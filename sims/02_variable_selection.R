# do variable selection

# variable selection based on all imputed training datasets --------------------
# @param dataset the stacked dataset consisting of all imputed datasets
# @param outcome_name the name of the outcome
# @param M the number of imputations
# @param k the cross-fitting iteration
# @param variable_selection_procedure the procedure (e.g., "lasso")
# @param impute_stability_threshold the proportion of imputed datasets in which
#   a variable must be selected to be part of the final set of variables
# @param alpha the desired per-comparison error rate
# @param sl_rank_threshold the rank-based threshold for SL variable selection;
#   by default,
# @param ss_b the number of SS iterations
# @param ss_pfer the SS parameter controlling the per-family error rate (expected num. false positives)
# @param ss_pi the SS parameter controlling the threshold (specify either this or \code{ss_q})
# @param ss_q the SS parameter controlling the number of variables considered (specify either this or \code{ss_pi})
# @param kf_fdr the desired knockoff FDR
# @param pfp_q the desired proportion of false positives (for intrinsic selection)
# @param gfwer_k the desired number of generalized family-wise errors we can make (for intrinsic selection)
# @param ss_b the number of bootstrap replications for stability selection
# @param ... other args (for Super Learner)
select_from_all <- function(imputed_training_data, outcome_name = "mucinous", M = 10, k = 1,
                            variable_selection_procedure = "lasso",
                            impute_stability_threshold = 0.7,
                            alpha = 0.04, sl_rank_threshold = 5, ss_pfer = 0.01,
                            ss_pi = 0.9, ss_q = 5, kf_fdr = 0.2, pfp_q = 0.75,
                            gfwer_k = 5, ss_b = 100, ...) {
  if (grepl("LJ", variable_selection_procedure) | grepl("lasso-BI-BL", variable_selection_procedure)) {
    # input is a list of bootstrap-imputed datasets
    if (grepl("LJ", variable_selection_procedure)) {
      # fit randomized lasso to each
      lambdas <- c(rev(seq(1, 5, by = 0.5)), rev(seq(0.025, 1, by = 0.1)))
      beta_list <- lapply(as.list(1:length(imputed_training_data)), function(b) {
        y <- imputed_training_data[[b]] %>% pull(!!outcome_name)
        x <- imputed_training_data[[b]] %>% select(-imp, -id, -!!outcome_name)
        weights <- runif(ncol(x), min = 0.5, max = 1)
        newx <- as.matrix(sweep(x, 2, weights, FUN = "/"))
        mod <- glmnet::glmnet(x = newx, y = y, lambda = lambdas, intercept = TRUE,
                              family = binomial())
        all_betas <- as.matrix(mod$beta)
        all_supports <- apply(all_betas, 2, function(col) as.numeric(abs(col) > 0))
        vim_select_tib <- data.table::rbindlist(lapply(as.list(1:ncol(all_betas)), function(j) {
          tibble::tibble(feature = 1:length(all_betas[, j]), importance = all_betas[, j],
                         lambda = mod$lambda[j], selected = all_supports[, j], b = b) 
        }))
      })
      all_betas <- data.table::rbindlist(beta_list)
      pooled_set <- all_betas %>%
        group_by(feature, lambda) %>%
        summarize(pi_lambda_j = mean(selected), .groups = "drop") %>%
        group_by(feature) %>%
        filter(pi_lambda_j == max(pi_lambda_j)) %>%
        filter(lambda == max(lambda)) %>%
        mutate(bi_ss_selected = as.numeric(pi_lambda_j >= ss_pi)) %>%
        pull(bi_ss_selected)
      all_selected_vars <- NULL
    } else {
      # fit CV lasso to each
      beta_list <- lapply(imputed_training_data, function(dat) {
        y <- dat %>% pull(!!outcome_name)
        x <- dat %>% select(-imp, -id, -!!outcome_name)
        feature_nms <- tibble(feature = names(x), feature_num = 1:ncol(x))
        mod <- glmnet::cv.glmnet(x = as.matrix(x), y = y, nfolds = 10, intercept = TRUE,
                               family = binomial())
        vim_select_tib <- flevr::extract_importance_glmnet(fit = mod, feature_names = names(x),
                                                          coef = 1) %>%
          left_join(feature_nms, by = "feature") %>%
          select(feature, importance, rank, feature_num) %>%
          rename(est = importance) %>%
          arrange(feature_num) %>%
          mutate(selected = as.numeric(abs(est) >= 0))
      })
      # obtain final estimated active set
      all_selected_vars <- do.call(rbind, lapply(
        beta_list, function(beta) {
          beta$selected
        }
      ))
      pooled_set <- as.numeric(colMeans(all_selected_vars) >= ss_pi)
    }
    # estimate betas from selected set
    est_vim <- data.table::rbindlist(lapply(as.list(1:length(imputed_training_data)), function(b) {
      y <- imputed_training_data[[b]] %>% pull(!!outcome_name)
      x <- imputed_training_data[[b]] %>% select(-imp, -id, -!!outcome_name)
      coefs <- rep(0, ncol(x))
      newx <- as.matrix(x[, pooled_set == 1])
      mod <- glm(y ~ newx, family = binomial())
      coefs[pooled_set == 1] <- coef(mod)[-1]
      data.frame(feature = 1:ncol(x), coef = coefs)
    })) %>%
      group_by(feature) %>%
      summarize(point_est = mean(coef), var = var(coef)) %>%
      mutate(est = ifelse(var > 0, abs(point_est) / var, 0))
  } else if (grepl("RR", variable_selection_procedure)) {
    # estimate SPVIMs for each imputed dataset
    spvim_lst <- lapply(as.list(seq_len(M)), function(m) {
      these_imputed_training_data <- imputed_training_data %>% 
        filter(imp == m) %>% 
        select(-imp, -id)
      y <- these_imputed_training_data %>% pull(!!outcome_name)
      x <- these_imputed_training_data %>% select(-!!outcome_name)
      est_spvim(x = x, y = y, n = nrow(x), ...)
    })
    control_lists <- list(
      list(quantity = "gFWER", base_method = "Holm", k = gfwer_k),
      list(quantity = "PFP", base_method = "Holm", k = gfwer_k, q = pfp_q),
      list(quantity = "FDR", base_method = "Holm", k = gfwer_k, q = pfp_q)
    )
    one_imputed_training_x <- imputed_training_data %>% 
      filter(imp == 1) %>% 
      select(-imp, -id, -!!outcome_name)
    # pool results based on Rubin's Rules, do variable selection based on the 
    # pooled VIMs
    selected_lst <- lapply(control_lists, function(control) {
      flevr::intrinsic_selection(spvim_ests = spvim_lst, sample_size = nrow(one_imputed_training_x),
                                 alpha = alpha, feature_names = names(one_imputed_training_x),
                                 control = control)
    })
    pooled_set <- lapply(selected_lst, function(l) l %>% pull(selected))
    est_vim <- selected_lst[[1]] %>%
      mutate(v = k)
    all_selected_vars <- NULL
  } else {
    # do variable selection for each imputed dataset separately
    selected_vars_lst <- vector("list", length = M)
    for (m in seq_len(M)) {
      these_imputed_training_data <- imputed_training_data %>%
        filter(imp == m) %>%
        select(-imp, -id)
      selected_vars_lst[[m]] <- select_from_one(dataset = these_imputed_training_data,
                                                outcome = outcome_name,
                                                variable_selection_procedure = variable_selection_procedure,
                                                alpha = alpha, sl_rank_threshold = sl_rank_threshold,
                                                ss_pfer = ss_pfer, ss_pi = ss_pi, ss_q = ss_q,
                                                kf_fdr = kf_fdr, pfp_q = pfp_q, gfwer_k = gfwer_k,
                                                ss_b = ss_b, ...)
    }
    all_selected_vars <- lapply(selected_vars_lst, function(l) l$selected_set)
    # estimated variable importance
    est_vim <- tibble::as_tibble(data.table::rbindlist(
      lapply(selected_vars_lst, function(l) l$est_vim)
    )) %>%
      mutate(m = rep(1:M, each = ncol(imputed_training_data) - 3),
             v = k)
    # pool the selected sets together
    if (grepl("SPVIM", variable_selection_procedure)) {
      gfwer_sets <- lapply(all_selected_vars, function(l) l[[1]])
      pfp_sets <- lapply(all_selected_vars, function(l) l[[2]])
      fdr_sets <- lapply(all_selected_vars, function(l) l[[3]])
      pooled_set <- list(
        flevr::pool_selected_sets(sets = gfwer_sets, threshold = impute_stability_threshold),
        flevr::pool_selected_sets(sets = pfp_sets, threshold = impute_stability_threshold),
        flevr::pool_selected_sets(sets = fdr_sets, threshold = impute_stability_threshold)
      )
    } else {
      pooled_set <- flevr::pool_selected_sets(sets = all_selected_vars,
                                              threshold = impute_stability_threshold)
    }
  }
  return(list(pooled = pooled_set, all = all_selected_vars,
              vim = est_vim))
}
# run variable selection for a single imputed dataset --------------------------
# @param dataset the dataset
# @param outcome the outcome of interest
# @param variable_selection_procedure the procedure (e.g., "lasso")
# @param alpha the desired per-comparison error rate
# @param sl_rank_threshold the rank-based threshold for SL variable selection;
#   by default,
# @param ss_b the number of SS iterations
# @param ss_pfer the SS parameter controlling the per-family error rate (expected num. false positives)
# @param ss_pi the SS parameter controlling the threshold (specify either this or \code{ss_q})
# @param ss_q the SS parameter controlling the number of variables considered (specify either this or \code{ss_pi})
# @param kf_fdr the desired knockoff FDR
# @param pfp_q the desired proportion of false positives (for intrinsic selection)
# @param gfwer_k the desired number of generalized family-wise errors we can make (for intrinsic selection)
select_from_one <- function(dataset, outcome, variable_selection_procedure = "lasso",
                                  alpha = 0.04, sl_rank_threshold = 5, ss_pfer = 0.01,
                                  ss_pi = 0.9, ss_q = 5, kf_fdr = 0.2, pfp_q = 0.75,
                                  gfwer_k = 5, ...) {
  y <- dataset %>% pull(!!outcome)
  x <- dataset %>% select(-!!outcome)
  feature_nms <- tibble(feature = names(x), feature_num = 1:ncol(x))
  if (grepl("lasso", variable_selection_procedure)) {
    if (!grepl("SS", variable_selection_procedure) & !grepl("KF", variable_selection_procedure)) {
      # run vanilla CV lasso
      mod <- glmnet::cv.glmnet(x = as.matrix(x), y = y, nfolds = 10, intercept = TRUE,
                               family = binomial())
      vim_select_tib <- flevr::extract_importance_glmnet(fit = mod, feature_names = names(x),
                                                         coef = 1) %>%
        left_join(feature_nms, by = "feature") %>%
        select(feature, importance, rank, feature_num) %>%
        rename(est = importance) %>%
        arrange(feature_num)
      est_vim <- vim_select_tib
      selected_set <- as.numeric(vim_select_tib$est != 0)
    } else if (grepl("SS", variable_selection_procedure)) {
      # run SS
      folds <- ss_folds(y = y, x = x, b = ss_b)
      mod <- stabs::stabsel(x = x, y = y, fitfun = stabs::glmnet.lasso,
                            args.fitfun = list(family = binomial(), type = "conservative"),
                            assumption = "none", verbose = FALSE, eval = TRUE, papply = mclapply,
                            mc.preschedule = FALSE, mc.cores = 1, q = ss_q, PFER = ss_pfer,
                            B = ss_b / 2, sampling.type = "SS", folds = folds)
      if (mod$cutoff > 0.9) {
        mod <- stabs::stabsel(mod, cutoff = 0.9)
      }
      est_vim <- tibble::tibble(
        feature = names(x), est = mod$max, rank = rank(-est), feature_num = 1:ncol(x)
      )
      selected_set <- rep(0, ncol(x))
      selected_set[mod$selected] <- 1
    } else {
      # run knockoffs
      x_mat <- as.matrix(x)
      check <- (glmnet:::weighted_mean_sd(x_mat)$sd == 0)
      if (any(check)) {
          x_mat <- x_mat[, -which(check)]
      }
      second_order_knockoffs <- knockoff::create.second_order(x_mat, method = "asdp", shrink = FALSE)
      kf_lasso_W <- knockoff::stat.glmnet_coefdiff(X = x_mat, X_k = second_order_knockoffs,
                                                   y = y, family = binomial(), cores = 1)
      kf_lasso_thresh <- knockoff::knockoff.threshold(kf_lasso_W, fdr = kf_fdr, offset = 1)
      if (is.infinite(kf_lasso_thresh)) {
        kf_lasso_thresh <- knockoff::knockoff.threshold(kf_lasso_W, fdr = kf_fdr, offset = 0)
      }
      if (any(check)) {
          tmp <- vector("numeric", length = ncol(x))
          tmp[!check] <- kf_lasso_W
          tmp[check] <- 0
          kf_lasso_W <- tmp
      }
      est_vim <- tibble::tibble(
        feature = names(x), est = kf_lasso_W, rank = rank(-abs(kf_lasso_W)),
        feature_num = 1:ncol(x)
      )
      selected_set <- as.numeric(kf_lasso_W >= kf_lasso_thresh)
    }
  } else if (grepl("SL", variable_selection_procedure)) {
    if (grepl("base", variable_selection_procedure)) {
      # no selection
      mod <- NULL
      est_vim <- tibble::tibble(feature = paste0("V", 1:p),
                                est = NA, rank = NA, feature_num = 1:p)
      selected_set <- rep(1, p)
    } else {
      args <- list(...)
      args$univariate_SL.library <- NULL
      if (!grepl("SS", variable_selection_procedure)) {
        # run SL
        mod <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = args$SL.library,
                                          family = args$family, method = args$method,
                                          cvControl = args$cvControl)
        vim_select_tib <- flevr::extrinsic_selection(fit = mod,
                                                     feature_names = names(x),
                                                     threshold = sl_rank_threshold,
                                                     import_type = "all", x = x, y = y) %>%
          left_join(feature_nms, by = "feature") %>%
          mutate(est = rank) %>%
          select(feature, est, rank, feature_num, selected) %>%
          arrange(feature_num)
        selected_set <- as.numeric(vim_select_tib$selected)
        est_vim <- vim_select_tib %>% select(-selected)
      } else {
        # run SL + SS
        folds <- ss_folds(y = y, x = x, b = ss_b)
        mod <- stabs::stabsel(x = x, y = y, fitfun = flevr::SL_stabs_fitfun,
                              args.fitfun = list(
                                SL.library = args$SL.library, family = args$family,
                                method = args$method, cvControl = args$cvControl
                              ),
                              assumption = "none", verbose = FALSE, eval = TRUE, papply = mclapply,
                              mc.preschedule = FALSE, mc.cores = 1, q = ss_q, PFER = ss_pfer,
                              B = ss_b / 2, sampling.type = "SS", folds = folds)
        if (mod$cutoff > 0.9) {
          mod <- stabs::stabsel(mod, cutoff = 0.9)
        }
        est_vim <- tibble::tibble(
          feature = names(x), est = mod$max, rank = rank(-est), feature_num = 1:ncol(x)
        )
        selected_set <- rep(0, ncol(x))
        selected_set[mod$selected] <- 1
      }
    }
  } else if (grepl("SPVIM", variable_selection_procedure)) {
    # set up
    control_lists <- list(
      list(quantity = "gFWER", base_method = "Holm", k = gfwer_k),
      list(quantity = "PFP", base_method = "Holm", k = gfwer_k, q = pfp_q),
      list(quantity = "FDR", base_method = "Holm", k = gfwer_k, q = pfp_q)
    )
    est_spvims <- est_spvim(x = x, y = y, n = nrow(dataset), ...)
    mod <- est_spvims$preds_lst
    est_vim <- est_spvims$mat %>%
      mutate(feature_num = as.numeric(s), rank = rank(-abs(est))) %>%
      select(est, p_value, rank, feature_num) %>%
      left_join(feature_nms, by = "feature_num") %>%
      select(feature, est, rank, p_value, feature_num)
    selected_lst <- lapply(control_lists, function(control) {
      flevr::intrinsic_selection(spvim_ests = est_spvims, sample_size = nrow(x),
                                 alpha = alpha, feature_names = names(x),
                                 control = control)
    })
    selected_set <- lapply(selected_lst, function(l) as.numeric(l$selected))
  }
  list(selected_set = selected_set, est_vim = est_vim)
}
