# generally useful functions

# Rubin's rules ----------------------------------------------------------------
# pool all estimates
# @param performance_lst a list of estimated performance
# @param variable_selection_procedure the variable selection procedure
# @param k the CV fold
# @param n the sample size
pool_all_ests <- function(performance_lst, variable_selection_procedure, k,
                          n) {
  if (grepl("SPVIM", variable_selection_procedure)) {
    pooled_lst <- lapply(as.list(c("gFWER", "PFP", "FDR")), function(proc) {
      all_ests <- unlist(lapply(performance_lst, function(l) l %>% filter(type == proc) %>%
                                  pull(perf)))
      all_vars <- unlist(lapply(performance_lst, function(l) l %>% filter(type == proc) %>%
                                  pull(var)))
      pool_ests(all_ests, all_vars)
    })
    pooled_tib <- as_tibble(rbindlist(pooled_lst))
    perf_tib <- tibble::tibble(id = k, est = pooled_tib$est,
                                    var = pooled_tib$total_var,
                                    n = n,
                                    type = c("gFWER", "PFP", "FDR"))
  } else {
    all_ests <- unlist(lapply(performance_lst, function(l) l$perf))
    all_vars <- unlist(lapply(performance_lst, function(l) l$var))
    pooled_lst <- pool_ests(all_ests, all_vars)
    perf_tib <- tibble::tibble(id = k, est = pooled_lst$est, var = pooled_lst$total_var,
                                    n = n)
  }
  return(perf_tib)
}
# pool estimates, ses using Rubin's rules
# @param ests a vector of estimates from M imputed datasets
# @param vars a vector of variance estimates from M imputed datasets
pool_ests <- function(ests, vars) {
  M <- length(ests)
  est <- mean(ests, na.rm = TRUE)
  U <- mean(vars, na.rm = TRUE)
  B <- M / (M - 1) * mean((ests - est)^2, na.rm = TRUE)
  total_var <- U + (1 + 1 / M) * B
  list(est = est, U = U, B = B, total_var = total_var)
}

# Create bootstrap datasets for Long + Johnson ---------------------------------
# @param data the original dataset
# @param B the number of bootstrap datasets
# @return a list of length B with the bootstrap datasets
generate_bootstrapped_data <- function(data, B = 100) {
  boot_data <- lapply(as.list(seq_len(B)), function(i) {
    data[sample(seq_len(nrow(data)), nrow(data), replace = TRUE), ]
  })
  return(boot_data)
}

# trim columns with more than desired proportion of missing data
# @param datset the dataset of interest
# @param missing_prop the proportion of missing data to trim if >
trim_na_cols <- function(dataset, missing_prop = 0.2) {
  if (missing_prop > 0) {
    dataset %>%
      purrr::discard(~ sum(is.na(.x)) / length(.x) * 100 > missing_prop)
  } else {
    dataset
  }
}

# trim columns with high pairwise correlation
trim_high_pairwise_corr <- function(dataset, corr = 0.9, outcome_name = "mucinous") {
  if (corr > 0) {
    tib_for_corr <- dataset %>% select(-any_of(c("id", "institution", outcome_name)))
    columns_to_remove <- tib_for_corr %>%
      as.matrix() %>%
      cor() %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      pivot_longer(-rowname, names_to = "colname", values_to = "correlation") %>%
      mutate(var_order = paste(rowname, colname) %>%
               strsplit(split = ' ') %>%
               map_chr( ~ sort(.x) %>%
                          paste(collapse = ' '))) %>%
      mutate(cnt = 1) %>%
      group_by(var_order) %>%
      mutate(cumsum = cumsum(cnt)) %>%
      filter(cumsum != 2) %>%
      ungroup() %>%
      select(-var_order, -cnt, -cumsum) %>%
      filter(!(rowname == colname)) %>%
      filter(abs(correlation) > corr)

    if (nrow(columns_to_remove) > 0) {
      dataset %>%
        select(-any_of(columns_to_remove$colname))
    } else {
      dataset
    }
  } else {
    dataset
  }
}

# find binary call thresholds
# @param score the score of interest
# @param call the binary call corresponding to the score
find_binary_call <- function(score, call, minimum = TRUE) {
  if (minimum) {
    min(score[call > 0], na.rm = TRUE)
  } else {
    max(score[call > 0], na.rm = TRUE)
  }
}

# prediction performance -------------------------------------------------------
# get CV-AUC
# @param preds the predicted values
# @param y the outcome values
# @param missing_threshold the threshold at which to return CV-AUC of non-NA observations or NA
get_cvauc <- function(preds, y, missing_threshold = 0.1) {
  perc_miss <- mean(is.na(preds))
  if (perc_miss <= missing_threshold) {
    cvauc <- tryCatch(
        cvAUC::ci.cvAUC(predictions = preds[!is.na(preds)], labels = y[!is.na(preds)]),
        error = function(e) list(cvAUC = NA, se = NA, ci = c(NA, NA), confidence = NA)
    )
  } else {
    cvauc <- list(cvAUC = NA, se = NA, ci = c(NA, NA), confidence = NA)
  }
  cvauc
}

# get SL preds
# @param newdata the new data
# @param fit the SL fit
get_sl_preds <- function(newdata, fit) {
  preds <- rep(NA, length = nrow(newdata))
  any_nas <- (rowSums(is.na(newdata)) > 0)
  if (any(any_nas)) {
    preds[!any_nas] <- predict(fit, newdata = newdata %>% filter(!any_nas), onlySL = TRUE)$pred
  } else {
    preds <- predict(fit, newdata = newdata, onlySL = TRUE)$pred
  }
  preds
}

# estimate prediction performance in a subgroup using cross-fitting
# @param imputed_datasets the imputed datasets
# @param M the number of imputations
# @param K the number of cross-fitting folds
# @param outcome_name the outcome
# @param variable_sets the selected variables
# @param subgroup the subgroup variable of interest
# @param n the sample size of the given subgroup
# @param prediction_procedure the prediction procedure
# @param arg_lst args to pass to prediction procedure
# @param missing_threshold threshold of remaining missing data to remove
subgroup_perf <- function(imputed_datasets, M = 10, K = 5, outcome_name = "mucinous",
                          variable_sets, subgroup = "overall", n,
                          selection_type,
                          prediction_procedure = "glm", arg_lst,
                          missing_threshold = 0.1) {
  # set up cross-fitting folds
  if (subgroup == "overall") {
    cross_fitting_folds <- vimp::make_folds(y = imputed_datasets %>% filter(imp == 1) %>% pull(outcome_name),
                                            V = K, stratified = TRUE)
  } else {
    one_imp <- imputed_datasets %>%
      filter(imp == 1)
    y <- one_imp %>% pull(!!outcome_name)
    cross_fitting_folds_0 <- vimp::make_folds(y = y[one_imp[, subgroup] == 0],
                                              V = K, stratified = TRUE)
    cross_fitting_folds_1 <- vimp::make_folds(y = y[one_imp[, subgroup] == 1],
                                              V = K, stratified = TRUE)
    cross_fitting_folds <- vector("numeric", n)
    cross_fitting_folds[one_imp[, subgroup] == 0] <- cross_fitting_folds_0
    cross_fitting_folds[one_imp[, subgroup] == 1] <- cross_fitting_folds_1
  }
  perf_lst <- vector("list", K)
  for (k in seq_len(K)) {
    if (subgroup == "overall") {
      train <- imputed_datasets %>%
        group_by(imp) %>%
        filter(cross_fitting_folds != k) %>%
        ungroup()
      test <- imputed_datasets %>%
        group_by(imp) %>%
        filter(cross_fitting_folds == k) %>%
        ungroup()
    } else {
      train <- imputed_datasets %>%
        group_by(imp) %>%
        filter(cross_fitting_folds != k) %>%
        ungroup() %>%
        rename(var = !!subgroup) %>%
        filter(var == 1)
      test <- imputed_datasets %>%
        group_by(imp) %>%
        filter(cross_fitting_folds == k) %>%
        ungroup() %>%
        rename(var = !!subgroup) %>%
        filter(var == 1)
      names(train) <- names(imputed_datasets)
      names(test) <- names(imputed_datasets)
    }
    if (n <= 1) {
      perf_lst[[k]] <- tibble::tibble(id = k, est = NA, var = NA, n = n,
                             type = NA, group = subgroup)
    } else {
      performance_lst <- prediction_performance(M = M, train = train,
                                                test = test, n = n,
                                                outcome_name = outcome_name,
                                                variable_sets = variable_sets,
                                                prediction_procedure = prediction_procedure,
                                                arg_lst = arg_lst, missing_threshold = missing_threshold)
      # combine using Rubin's rules
      perf_lst[[k]] <- pool_all_ests(performance_lst = performance_lst,
                            variable_selection_procedure = selection_type,
                            k = k, n = n) %>%
        mutate(group = subgroup)
    }
  }
  if (any(grepl("type", names(perf_lst[[1]])))) {
    perf <- as_tibble(rbindlist(perf_lst)) %>%
      group_by(n, type, group) %>%
      summarize(est = mean(est), var = mean(var), .groups = "drop")
  } else {
    perf <- as_tibble(rbindlist(perf_lst)) %>%
      group_by(n, group) %>%
      summarize(est = mean(est), var = mean(var), .groups = "drop") %>%
      mutate(type = NA)
  }
  return(perf)
}

# binary-outcome stability selection -------------------------------------------
# @param y the outcome
# @param x the covariates
# @param b the number of SS iterations
ss_folds <- function(y, x, b = 100) {
  if (length(unique(y)) == 2) {
    n_0 <- sum(y == 0)
    n_1 <- sum(y == 1)
    k_0 <- floor(n_0 * 0.5)
    k_1 <- ceiling(n_1 * 0.5)
    indices_0 <- rep(c(0, 1), c(n_0 - k_0, k_0))
    indices_1 <- rep(c(0, 1), c(n_1 - k_1, k_1))
    folds_0 <- replicate(b / 2, sample(indices_0))[sample(1:n_0), , drop = FALSE]
    folds_1 <- replicate(b / 2, sample(indices_1))[sample(1:n_1), , drop = FALSE]
    folds <- matrix(0, nrow = nrow(x), ncol = b / 2)
    folds[y == 0, ] <- folds_0
    folds[y == 1, ] <- folds_1
  } else {
    folds <- stabs::subsample(rep(1, length(y)), B = b / 2)
  }
  folds
}

# set up and estimate SPVIMs ---------------------------------------------------
# @param x the covariates
# @param y the outcome
# @param n the sample size
# @param gfwer_k the k for gFWER control
# @param pfp_q the q for PFP control
# @param ... other args to pass to Super Learner
est_spvim <- function(x, y, n, ...) {
  args <- list(...)
  this_v <- switch((n < 500) + 1, 2, 5)
  this_v_2 <- switch((n < 500) + 1, 2, 5)
  # estimate SPVIM
  gamma <- 1
  nrounds <- 10
  for (i in seq_len(nrounds)) {
    z_w_lst <- vimp::sample_subsets(p = ncol(x), n = nrow(x), gamma = gamma)
    if (nrow(z_w_lst$Z) > ncol(x) | (i == nrounds)) {
      break
    } else {
      gamma <- gamma + 0.5
    }
  }
  est_spvims <- vimp::sp_vim(Y = y, X = x, V = this_v, gamma = gamma,
                             type = "auc", stratified = TRUE,
                             SL.library = args$spvim_library,
                             univariate_SL.library = args$univariate_SL.library,
                             cvControl = list(V = this_v_2, stratifyCV = TRUE),
                             family = args$family, method = args$method)
  return(est_spvims)
}

# fitting Super Learners -------------------------------------------------------

method.CC_nloglik2 <- function () {
  computePred = function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      warning("All metalearner coefficients are zero, predictions will all be 0",
              call. = FALSE)
    }
    plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
             matrix(coef[coef != 0]))
  }
  computeCoef = function(Z, Y, libraryNames, obsWeights, control,
                         verbose, ...) {
    tol <- 8
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
                                                  ifelse(Y, plogis(x, log.p = TRUE), plogis(x, log.p = TRUE,
                                                                                            lower.tail = FALSE))))
    names(cvRisk) <- libraryNames
    obj_and_grad <- function(y, x, w = NULL) {
      y <- y
      x <- x
      function(beta) {
        xB <- x %*% cbind(beta)
        loglik <- y * plogis(xB, log.p = TRUE) + (1 -
                                                    y) * plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w))
          loglik <- loglik * w
        obj <- -2 * sum(loglik)
        p <- plogis(xB)
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
    r <- nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)),
                        eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds,
                        ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) -
                                                                         1), eval_jac_g_eq = function(beta) rep(1, length(beta)),
                        opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
    if (r$status < 1 || r$status > 4) {
      warning(r$message)
    }
    coef <- r$solution
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

# correlation screening
screen.corRank.p <- function(Y, X, family, method = "spearman",
                             rank = ifelse(ncol(X) > 2, ncol(X) /
                                             switch((ncol(X) > 100) + 1, 2, 10), 2),
                             minVar = 1, ...) {
  listp <- apply(X, 2, function(x, Y, method) {
    ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method, exact = FALSE)$p.value)
  }, Y = Y, method = method)
  whichVariable <- (rank(listp) <= rank)
  if (sum(whichVariable) < minVar) {
      the_variable <- sample(seq_len(ncol(X)), minVar)
      whichVariable[the_variable] <- TRUE
  }
  return(whichVariable)
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
                            eval_metric = "logloss",
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
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = TRUE, ...) {
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
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))){
  SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
  SL.ranger.imp(..., X  = X, mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = min(ncol(X), floor(sqrt(ncol(X)) * 2))){
  SL.ranger.imp(..., X = X, mtry = mtry)
}
# function used to do smarter CV for glmnet
get_fold_id <- function(Y){
  fold_id <- rep(0, length(Y))
  wiY0 <- which(Y == 0)
  wiY1 <- which(Y == 1)
  #if <4 cases, no cv
  if(length(wiY1) == 4){
    #if exactly 4 cases, 4-fold cv
    #1 case per fold
    fold <- 1:4
    fold_id[sample(wiY1)] <- fold
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }else{
    #if >=5 cases, 5 fold cv
    #cases split as evenly as possible
    fold <- 1:5
    fold_id[sample(wiY1)] <- rep(fold, length = length(wiY1))
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }
  return(fold_id)
}

# lasso
SL.glmnet.0 <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, alpha = 1, nfolds = 5,
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
  SuperLearner:::.SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }

  if(family$family == "binomial"){
    fold_id <- get_fold_id(Y)
    nfolds <- max(fold_id)
    if(nfolds != 0){
      fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                                 lambda = NULL, type.measure = loss, nfolds = nfolds,
                                 foldid = fold_id, family = family$family, alpha = alpha, nlambda = nlambda,
                                 ...)
      pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                        "lambda.min", "lambda.1se"))
      fit <- list(object = fitCV, useMin = useMin)
      class(fit) <- "SL.glmnet"
    }else{
      # if fewer than 3 cases, just use mean
      meanY <- weighted.mean(Y, w = obsWeights)
      pred <- rep.int(meanY, times = nrow(newX))
      fit <- list(object = meanY)
      out <- list(pred = pred, fit = fit)
      class(fit) <- c("SL.mean")
    }
  }else{
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                               lambda = NULL, type.measure = loss, nfolds = nfolds,
                               family = family$family, alpha = alpha, nlambda = nlambda,
                               ...)
    pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                      "lambda.min", "lambda.1se"))
    fit <- list(object = fitCV, useMin = useMin)
    class(fit) <- "SL.glmnet"
  }
  out <- list(pred = pred, fit = fit)
  return(out)
}
# elastic net
SL.glmnet.50 <- function(..., alpha = 0.5){
  SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.25 <- function(..., alpha = 0.25){
  SL.glmnet.0(..., alpha = alpha)
}
SL.glmnet.75 <- function(..., alpha = 0.75){
  SL.glmnet.0(..., alpha = alpha)
}
# polspline
my_SL.polymars <- function(Y, X, newX, family, obsWeights, ...) {
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
#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function ()
{
  computeCoef = function(Z, Y, libraryNames, verbose, obsWeights,
                         errorsInLibrary = NULL, ...) {
    cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x -
                                                           Y)^2))
    names(cvRisk) <- libraryNames
    compute <- function(x, y, wt = rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      fit <- tryCatch({quadprog::solve.QP(Dmat = D, dvec = d, Amat = A,
                                          bvec = bvec, meq = 1)
      }, error = function(e){
        out <- list()
        class(out) <- "error"
        out
      })
      invisible(fit)
    }
    modZ <- Z
    naCols <- which(apply(Z, 2, function(z) {
      all(z == 0)
    }))
    anyNACols <- length(naCols) > 0
    if (anyNACols) {
      warning(paste0(paste0(libraryNames[naCols], collapse = ", "),
                     " have NAs.", "Removing from super learner."))
    }
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if (anyDupCols) {
      warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                     " are duplicates of previous learners.", " Removing from super learner."))
    }
    if (anyDupCols | anyNACols) {
      rmCols <- unique(c(naCols, dupCols))
      modZ <- Z[, -rmCols, drop = FALSE]
    }
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    if(class(fit) != "error"){
      coef <- fit$solution
    }else{
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    if(class(fit) != "error"){
      if (anyDupCols | anyNACols) {
        ind <- c(seq_along(coef), rmCols - 0.5)
        coef <- c(coef, rep(0, length(rmCols)))
        coef <- coef[order(ind)]
      }
      coef[coef < 1e-04] <- 0
      coef <- coef/sum(coef)
    }
    if (!sum(coef) > 0)
      warning("All algorithms have zero weight", call. = FALSE)
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(require = "quadprog", computeCoef = computeCoef,
              computePred = computePred)
  invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
tmp_method.CC_nloglik <- function ()
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

# Summarize CV.SuperLearner objects and extract point estimates, CIs -----------
# this function takes as input a fitted object EITHER of class SuperLearner OR
# of class CV.SuperLearner and computes a summary table of performance.
# if all object$Y are 0/1, it will compute AUC; otherwise it computes R^2
# if SuperLearner is used to evaluate CV performance of a single algorithm,
# a single row table is returned with that algorithms performance.
# if CV.SuperLearner is used to evaluate CV performance of a CV-tuned single algorithm,
# a table with a row for each choice of tuning parameters and for the cv-selected tuning
# parameters is returned.
# if CV.SuperLearner is used to evaluated CV performance of SuperLearner, then an
# additional row is added that describes the performance of SuperLearner.
summary.myCV.SuperLearner <- function (object, obsWeights = NULL, method = NULL, opts, ...) {
  if ("env" %in% names(object)) {
    env = object$env
  }else {
    env = parent.frame()
  }

  is_sl <- "SuperLearner" %in% class(object)
  is_cvsl <- "CV.SuperLearner" %in% class(object)
  if(is_sl | is_cvsl){
    library.names <- object$libraryNames
    if(is_cvsl){
      V <- object$V
    }else{
      V <- length(object$validRows)
    }
    n <- length(object$SL.predict)
    if (is.null(obsWeights)) {
      obsWeights <- rep(1, length(object$Y))
    }

    if(is_cvsl){
      folds <- object$folds
    }else if(is_sl){
      folds <- object$validRows
    }

    if(is_cvsl){
      # only will use this if multiple learners selected
      SL.predict <- object$SL.predict
      # this will be only output if single learner used and opts$cvtune
      discreteSL.predict <- object$discreteSL.predict
      # only will use this if multiple learners selected
      library.predict <- object$library.predict
    }else if(is_sl){
      # in this case a single "default" learner was requested
      # so we can pull Z out from the object
      SL.predict <- object$Z[,1]
    }

    Y <- object$Y
    Risk.SL <- rep(NA, length = V)
    se.SL <- rep(NA, length = V)
    if(is_cvsl){
      Risk.dSL <- rep(NA, length = V)
      se.dSL <- rep(NA, length = V)
      Risk.library <- matrix(NA, nrow = length(library.names),
                             ncol = V)
      se.library <- matrix(NA, nrow = length(library.names),
                           ncol = V)
      rownames(Risk.library) <- library.names
    }
    if (!(all(Y %in% c(0,1)))) {
      for (ii in seq_len(V)) {
        Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                                                         SL.predict[folds[[ii]]])^2)
        if(is_cvsl){
          Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                                                            discreteSL.predict[folds[[ii]]])^2)
          Risk.library[, ii] <- apply(library.predict[folds[[ii]],
                                                      , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] *
                                                                                             (Y[folds[[ii]]] - x)^2))
        }
      }
      if_sl <- (Y - SL.predict)^2 - mean((Y - SL.predict)^2)
      if(is_cvsl){
        if_dsl <- (Y - discreteSL.predict)^2 - mean((Y - discreteSL.predict)^2)
        if_library <- apply(library.predict, 2, function(x){ (Y - x)^2 - mean((Y - x)^2) })
      }
      if_varY <- (Y - mean(Y))^2 - mean((Y - mean(Y))^2)
      get_log_se <- function(if_risk, if_varY, risk, varY,
                             n = length(if_risk)){
        grad <- matrix(c(1 / risk, - 1 /varY), nrow = 2)
        Sig <- cov(cbind(if_risk, if_varY))
        se_log <- t(grad) %*% Sig %*% grad
        return(se_log)
      }

      if(is_cvsl){
        if(length(opts$learners) > 1){
          se <- (1/sqrt(n)) * c(
            get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y)),
            get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
            mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                   function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
          )
          Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
                                            library.names), Ave = c(1 - mean(Risk.SL)/var(Y), 1 - mean(Risk.dSL)/var(Y),
                                                                    apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se, Min = c(min(1 - Risk.SL/var(Y)),
                                                                                                                                                     min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})), Max = c(max(1 - Risk.SL/var(Y)),
                                                                                                                                                                                                                                                       max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))

        }else{
          se <- (1/sqrt(n)) * c(
            get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
            mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                   function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
          )

          Table <- data.frame(Algorithm = c("Discrete SL",
                                            library.names), Ave = c(1 - mean(Risk.dSL)/var(Y),
                                                                    apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se,
                              Min = c(min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})),
                              Max = c(max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))
        }
      }else{
        se <- (1/sqrt(n)) * get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y))
        Table <- data.frame(Algorithm = c(library.names[1]), Ave = c(1 - mean(Risk.SL)/var(Y)),
                            log_se = se,
                            Min = c(min(1 - Risk.SL/var(Y))),
                            Max = c(max(1 - Risk.SL/var(Y))))
      }
    }else {
      requireNamespace("cvAUC")
      for (ii in seq_len(V)) {
        sl_auc <- cvAUC::ci.cvAUC(predictions = SL.predict[folds[[ii]]],
                                  labels = Y[folds[[ii]]], folds = NULL)
        Risk.SL[ii] <- sl_auc$cvAUC
        se.SL[ii] <- sl_auc$se
        if(is_cvsl){
          dsl_auc <- cvAUC::ci.cvAUC(predictions = discreteSL.predict[folds[[ii]]],
                                     labels = Y[folds[[ii]]], folds = NULL)
          Risk.dSL[ii] <- dsl_auc$cvAUC
          se.dSL[ii] <- dsl_auc$se
          library_auc <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x){
            tmp <- cvAUC::ci.cvAUC(predictions = x, labels = Y[folds[[ii]]], folds = NULL)
            return(c(tmp$cvAUC, tmp$se))
          })
          Risk.library[,ii] <- library_auc[1,]
          se.library[,ii] <- library_auc[2,]
        }
      }
      if(is_cvsl){
        if(length(opts$learners) > 1){
          se <- c(mean(se.SL, na.rm = TRUE), mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
          Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
                                            library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL),
                                                                    apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL),
                                                                                                                    min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL),
                                                                                                                                                                         max(Risk.dSL), apply(Risk.library, 1, max)))
        }else{
          se <- c(mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
          Table <- data.frame(Algorithm = c("Discrete SL",
                                            library.names), Ave = c(mean(Risk.dSL),
                                                                    apply(Risk.library, 1, mean)), se = se, Min = c(
                                                                      min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(
                                                                        max(Risk.dSL), apply(Risk.library, 1, max)))
        }
      }else{
        se <- c(mean(se.SL, na.rm = TRUE))
        Table <- data.frame(Algorithm = c(library.names[1]),
                            Ave = c(mean(Risk.SL)), se = se,
                            Min = c(min(Risk.SL)),
                            Max = c(max(Risk.SL)))
      }
    }
    out <- list(call = object$call, method = method, V = V, Table = Table)
  }
  class(out) <- "summary.myCV.SuperLearner"
  return(out)
}
