## utility functions

# take the expit
# @param x the argument to take the expit of
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

# make missing data patterns for mice::ampute
# @param type the type of missing data (e.g., "monotone" or "non-monotone")
# @param p the number of variables
make_missing_pattern <- function(type = "monotone", p = 50) {
  # set up missingness pattern:
  # 0 indicates a variable should have missing values, 1 complete
  if (type == "monotone") {
    main_vars_plus_outcome <- expand.grid(
      y = 1, x1 = 1, x2 = 0:1, x3 = 1, x4 = 0:1, x5 = 1, x6 = 0:1
    )[c(4, 2, 1), ]

  } else {
    main_vars_plus_outcome <- expand.grid(
      y = 1, x1 = 1, x2 = 0:1, x3 = 1, x4 = 0:1, x5 = 1, x6 = 0:1
    )[-8, ]
  }
  # missing pattern is for the first 6; observe most of the noise variables (approx 91%)
  if (p > 6) {
    missing_noise_vars <- matrix(0, nrow = nrow(main_vars_plus_outcome),
                                 ncol = ceiling((40 / (500 - 6)) * p))
    observed_noise_vars <- matrix(1, nrow = nrow(main_vars_plus_outcome),
                                  ncol = p - ncol(missing_noise_vars) - 6)

    missing_pattern <- cbind(main_vars_plus_outcome, observed_noise_vars,
                             missing_noise_vars)
  } else {
    missing_pattern <- main_vars_plus_outcome
  }
  missing_pattern
}

# boosted trees
SL.xgboost_new <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000,
    max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
    nthread = 1, verbose = 0, save_period = NULL, ...)
{
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror", ## if xgboost version >=1.1.1.1, changed from reg:linear to reg:squarederror
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
            nthread = nthread, params = params, save_period = save_period)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}

# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)),
                           num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = FALSE, ...) {
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose, importance = "impurity")
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))) {
  SL.ranger.imp(..., X = X, mtry = mtry)
}
SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)) {
  SL.ranger.imp(..., X = X, mtry = mtry)
}
SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)) {
  SL.ranger.imp(..., X = X, mtry = mtry)
}
# polspline
SL.polymars2 <- function(Y, X, newX, family, obsWeights, ...) {
  SuperLearner:::.SL.require("polspline")
  if (family$family == "gaussian") {
    capture.output(fit.mars <- polspline::polymars(Y, X, weights = obsWeights, verbose = FALSE))
    pred <- predict(fit.mars, x = newX)
    fit <- list(object = fit.mars)
  }
  if (family$family == "binomial") {
    capture.output(fit.mars <- polspline::polyclass(Y, X, cv = 5, delete = 0,
                                     weight = obsWeights, silent = TRUE,
                                     select = 1))
    pred <- polspline::ppolyclass(cov = newX, fit = fit.mars)[, 2]
    fit <- list(fit = fit.mars)
  }
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.polymars")
  return(out)
}
# glmnet screen
screen.glmnet_new <- function(Y, X, family, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100,  ...) {
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  # the [-1] removes the intercept
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    # only create a new cut if *any* of the lambdas will result in minscreen vars getting in
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- switch(as.numeric(any(sumCoef >= minscreen)) + 1,
                            (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0),
                            (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0))
  }
  if (sum(whichVariable) == 0) {
    whichVariable <- rep(TRUE, ncol(X))
  }
  return(whichVariable)
}
#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
method.CC_nloglik2 <- function ()
{
  computePred = function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      stop("All metalearner coefficients are zero, cannot compute prediction.")
    }
    stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
                    matrix(coef[coef != 0]))
  }
  computeCoef = function(Z, Y, libraryNames, obsWeights, control,
                         verbose, ...) {
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    modZ <- Z
    if (anyDupCols) {
      warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                     " are duplicates of previous learners.", " Removing from super learner."))
      modZ <- modZ[, -dupCols, drop = FALSE]
    }
    modlogitZ <- trimLogit(modZ, control$trimLogit)
    logitZ <- trimLogit(Z, control$trimLogit)
    cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights *
                                                  ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x, log.p = TRUE,
                                                                                                          lower.tail = FALSE))))
    names(cvRisk) <- libraryNames
    obj_and_grad <- function(y, x, w = NULL) {
      y <- y
      x <- x
      function(beta) {
        xB <- x %*% cbind(beta)
        loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
                                                           y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w))
          loglik <- loglik * w
        obj <- -2 * sum(loglik)
        p <- stats::plogis(xB)
        grad <- if (is.null(w))
          2 * crossprod(x, cbind(p - y))
        else 2 * crossprod(x, w * cbind(p - y))
        list(objective = obj, gradient = grad)
      }
    }
    lower_bounds = rep(0, ncol(modZ))
    upper_bounds = rep(1, ncol(modZ))
    if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] = 0
    }
    r <- tryCatch({nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)),
                                  eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds,
                                  ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) -
                                                                                   1), eval_jac_g_eq = function(beta) rep(1, length(beta)),
                                  opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
    }, error = function(e){
      out <- list()
      class(out) <- "error"
      out
    })
    if (r$status < 1 || r$status > 4) {
      warning(r$message)
    }
    if(class(r) != "error"){
      coef <- r$solution
    }else{
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (anyDupCols) {
      ind <- c(seq_along(coef), dupCols - 0.5)
      coef <- c(coef, rep(0, length(dupCols)))
      coef <- coef[order(ind)]
    }
    coef[coef < 1e-04] <- 0
    coef <- coef/sum(coef)
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}

## utility functions for super learner
get_sl_opts <- function(family) {
  if (family == "gaussian") {
    sl_fam <- gaussian()
    sl_method <- "method.CC_LS"
  } else {
    if (grepl("probit", family)) {
      sl_fam <- binomial(link = "probit")
    } else {
      sl_fam <- binomial()
    }
    sl_method <- "method.CC_nloglik2"
  }
  return(list(fam = sl_fam, method = sl_method))
}

# get the procedure, base method, and control quantity for SPVIM/SAGE+-based variable selection
get_vim_proc <- function(est_type) {
  rem_str <- paste0(switch((grepl("SPVIM", est_type)) + 1, "SAGE+", "SPVIM"), "-")
  # get the procedure and base method
  proc <- gsub(rem_str, "", est_type)
  split_proc <- unlist(strsplit(proc, "-", fixed = TRUE))
  return(list(procedure = split_proc[1], quantity = split_proc[1],
              base_method = split_proc[2], fdr_method = split_proc[2]))
}

# get the selected set and set size after return (since it returns as a string, for ease)
get_selected_sets <- function(x) {
  split_str <- stringr::str_split(x, ";")
  split_num <- lapply(split_str, as.numeric)
  return(split_num)
}
get_selected_set_size <- function(x) {
  selected_sets <- get_selected_sets(x)
  return(unlist(lapply(selected_sets, length)))
}

# get the "coefs" (or equivalent) for a fit object
#' @param fit the fitted object
#' @return a character vector
get_fit_out <- function(fit) {
  if (any(grepl("glmnet", class(fit)))) {
    coefs <- coef(fit$glmnet.fit, s = fit$lambda.min)[-1]
    fit_out <- paste0(names(coefs), ": ", round(coefs, 3), collapse = "; ")
  } else if (any(grepl("glm", class(fit)))) {
    coefs <- fit$coefficients
    fit_out <- paste0(names(coefs), ": ", round(coefs, 3), collapse = "; ")
  } else if (any(grepl("SuperLearner", class(fit)))){
    out_df <- rbind.data.frame(round(fit$cvRisk, 3), round(fit$coef, 3))
    names(out_df) <- names(fit$coef)
    fit_out <- paste0(
      names(out_df), ": ", apply(out_df, 2, function(x) paste0(x, collapse = "/")),
      collapse = "; "
    )
  } else {
    fit_out <- ""
  }
  fit_out
}
# get the performance of the selected set
#' @param selected_set a vector of 0s (not selected) and 1s (selected)
#'                     corresponding to the number of predictors
#' @param x the covariates
#' @param test_x the test-set covariates
#' @param test_y the test-set outcome
#' @return
get_selected_set_perf <- function(selected_set = NULL, x = NULL, y = NULL,
                                  test_x = NULL, test_y = NULL,
                                  estimator_type = "lasso",
                                  learner_lib = c("SL.ranger", "SL.mean"),
                                  family = "binomial") {
  if (sum(selected_set) == 1) {
    # add a column of zeros to trick glmnet
    selected_vars <- as.matrix(cbind(rnorm(nrow(x), 0, 1e-4), x[, selected_set == 1, drop = FALSE]))
    colnames(selected_vars) <- paste0("V", 0:1)
    test_selected_vars <- as.matrix(cbind(0, test_x[, selected_set == 1,
                                                    drop = FALSE]))
    colnames(test_selected_vars) <- paste0("V", 0:1)
  } else {
    selected_vars <- as.matrix(x[, selected_set == 1, drop = FALSE])
    test_selected_vars <- as.matrix(test_x[, selected_set == 1, drop = FALSE])
  }
  if (sum(selected_set) == 0) {
    selected_mod <- NA
    selected_mod_preds <- rep(mean(test_y), length(test_y))
  } else {
    opts <- get_sl_opts(family)
    if (grepl("lasso", estimator_type)) {
      # selected_mod <- glmnet::cv.glmnet(x = selected_vars, y = y,
      #                                   nfolds = 10, intercept = TRUE,
      #                                   family = opts$fam)
      # selected_mod_preds <- predict(selected_mod, s = "lambda.min",
      #                               newx = test_selected_vars)
      selected_mod <- glm(y ~ ., data = data.frame(y = y, selected_vars),
                          family = opts$fam)
      selected_mod_preds <- predict(selected_mod, newdata = data.frame(test_selected_vars),
                                    type = switch((opts$fam$family == "binomial") + 1, "link", "response"))
    } else {
      selected_mod <- SuperLearner::SuperLearner(
        Y = y, X = as.data.frame(selected_vars), SL.library = learner_lib,
        cvControl = list(V = 5,
                         stratifyCV = (opts$fam$family == "binomial")), family = opts$fam,
        method = opts$method
      )
      selected_mod_preds <- predict(
        selected_mod, newdata = as.data.frame(test_selected_vars), onlySL = TRUE
      )$pred
    }
  }
  if (grepl("binomial", family)) {
    selected_mod_perf <- cvAUC::cvAUC(selected_mod_preds, test_y)$cvAUC
  } else {
    selected_mod_perf <- vimp::measure_r_squared(selected_mod_preds,
                                                 test_y)$point_est
  }
  fit_out <- get_fit_out(selected_mod)
  list(perf = selected_mod_perf, fit_out = fit_out)
}

## get sensitivity
#' @param selected_sets vector of selected sets (characters)
#' @param true_active_set the true active set
#' @return sensitivity (true positive rate, i.e., proportion of active set recovered)
get_sensitivity <- function(selected_sets, true_active_set) {
  split_str <- selected_sets %>%
    stringr::str_split(pattern = ";")
  all_set_lst <- lapply(split_str, as.numeric)
  s0 <- length(true_active_set)
  sens <- unlist(lapply(all_set_lst,
                        function(x) length(dplyr::intersect(true_active_set, x)))) / s0
  return(sens)
}

# set specificity for SPVIM based on n
#' @param n the sample size
#' @param ns the sample sizes
#' @param p the dimension
#' @param max_spec the maximum specificity
#' @return the target specificity
set_specificity <- function(n, ns = c(200, 500, 1000, 2000, 3000), p = 30, max_spec = 0.9) {
  # target at n = 200 is 0.7, at n = 3000 is 0.9
  if (p == 30) {
    spec <- seq(max_spec - 0.1, max_spec, length.out = length(ns))
  } else {
    spec <- seq(max_spec - 0.1, max_spec, length.out = length(ns))
  }
  fit <- suppressWarnings(glm(spec ~ ns, family = binomial))
  min(predict(fit, newdata = data.frame(ns = n), type = "response"), 1)
}

## get specificity
#' @param selected_sets vector of selected sets (characters)
#' @param true_active_set the true active set
#' @param p the dimension
#' @return specificity (true negative rate, i.e., proportion of non-active set not recovered)
get_specificity <- function(selected_sets, true_active_set, p) {
  split_str <- selected_sets %>%
    stringr::str_split(pattern = ";")
  all_set_lst <- lapply(split_str, as.numeric)
  true_nonactive_sets <- lapply(as.list(p), function(x) dplyr::setdiff(1:x, true_active_set))
  s0 <- lapply(true_nonactive_sets, length)
  spec <- mapply(function(x, y, z) length(dplyr::setdiff(x, y))/z,
                 x = true_nonactive_sets, y = all_set_lst, z = s0)
  return(spec)
}

# obtain q for stability selection from PFER_max and pi_thr
#' @param PFER_max the maximum per-family error rate (generally p / 2 * per-comparison error rate (defaults to 0.05))
#' @param pi_thr the threshold value for stability
#' @return q the number of features selected by SS per boosting run
get_ss_q <- function(PFER_max = 0.375, pi_thr = 0.75, p = 15) {
  sqrt(p) * sqrt((2 * pi_thr - 1) * PFER_max)
}

.make_folds <- function(y, V, stratified = FALSE, probs = rep(1/V, V)) {
    folds <- vector("numeric", length(y))
    if (length(unique(folds)) == 1) {
      if (stratified) {
          folds_1 <- sample(rep(seq_len(V), length = sum(y == 1)))
          folds_0 <- sample(rep(seq_len(V), length = sum(y == 0)))
          folds[y == 1] <- folds_1
          folds[y == 0] <- folds_0
      } else {
          folds <- sample(rep(seq_len(V), length = length(y)))
      }
    } else {
      if (stratified) {
        folds_1 <- rep(seq_len(V), probs * sum(y == 1))
        folds_1 <- c(folds_1, sample(seq_len(V), size = sum(y == 1) - length(folds_1),
                                     replace = TRUE, prob = probs))
        folds_0 <- rep(seq_len(V), probs * sum(y == 0))
        folds_0 <- c(folds_0, sample(seq_len(V), size = sum(y == 0) - length(folds_0),
                                     replace = TRUE, prob = probs))
        folds_1 <- sample(folds_1)
        folds_0 <- sample(folds_0)
        folds[y == 1] <- folds_1
        folds[y == 0] <- folds_0
      } else {
        folds <- rep(seq_len(V), probs * length(y))
        folds <- c(folds, sample(seq_len(V), size = length(y) - length(folds),
                                 replace = TRUE, prob = probs))
        folds <- sample(folds)
      }
    }
    return(folds)
}

# screen based on rank correlation
screen.corRank.p <- function(Y, X, family, method = "spearman",
                             rank = ifelse(ncol(X) > 2, ncol(X) /
                               switch((ncol(X) > 100) + 1, 2, 10), 2),
                             minVar = 1, ...) {
  listp <- apply(X, 2, function(x, Y, method) {
    ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method, exact = FALSE)$p.value)
  }, Y = Y, method = method)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}
screen.corRank.10 <- function(Y, X, family, method = "spearman", rank = 10, ...) {
  screen.corRank(Y = Y, X = X, family = family, method = method, rank = rank, ...)
}
screen.corRank.25 <- function(Y, X, family, method = "spearman", rank = 25, ...) {
  screen.corRank(Y = Y, X = X, family = family, method = method, rank = rank, ...)
}
screen.corRank.50 <- function(Y, X, family, method = "spearman", rank = 50, ...) {
  screen.corRank(Y = Y, X = X, family = family, method = method, rank = rank, ...)
}

# make estimator levels
make_est_levels <- function(tib) {
  levels_init <- unique(paste(tib$selection_type, tib$extra_layer, sep = " + "))
  remove_none <- gsub(" + none", "", gsub("_none + none", "", levels_init, fixed = TRUE), fixed = TRUE)
  final_levels <- gsub("_BH", " + rank (BH)", gsub("_none + screen", " (screen)",
                                                   gsub("_screen", " + screened rank",
                                                        gsub("_rank", " + rank",
                                                             gsub("sl", "SL",
                                                                  gsub("spvim", "SPVIM",
                                                                       gsub("sage", "SAGE+",
                                                                            remove_none, fixed = TRUE))))),
                                                   fixed = TRUE), fixed = TRUE)
  final_levels
}
