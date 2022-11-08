# estimate prediction performance based on the chosen algorithm

# get prediction performance based on all sets of variables
# @param M the number of imputations
# @param train the stacked dataset of training data including all imputations
# @param test the stacked dataset of testing data including all imputations
# @param variable_sets a list with the correct set of variables to use for each imputed dataset
# @param n the overall sample size
# @param prediction_procedure the prediction procedure to use (glm, SL)
# @param arg_lst a list of arguments to pass to Super Learner
# @param missing_threshold the missing data threshold
prediction_performance <- function(M = 10, train, test, variable_sets, n,
                                   outcome_name = "mucinous", prediction_procedure = "glm",
                                   arg_lst, missing_threshold) {
    performance_lst <- vector("list", length = M)
    for (m in seq_len(M)) {
        these_selected_vars <- variable_sets[[m]]
        # get train/test for *this* imputation
        these_imputed_training_data <- train %>%
          filter(imp == m) %>%
          select(-imp, -id)
        these_imputed_testing_data <- test %>%
          filter(imp == m) %>%
          select(-imp, -id)
        # split into outcome and covariates
        y_train <- these_imputed_training_data %>% pull(!!outcome_name)
        x_train <- these_imputed_training_data %>% select(-!!outcome_name)
        y_test <- these_imputed_testing_data %>% pull(!!outcome_name)
        x_test <- these_imputed_testing_data %>% select(-!!outcome_name)
        # train and assess performance on test
        if (is.list(these_selected_vars)) {
            cvaucs <- lapply(these_selected_vars, function(vars) {
                one_pred_perf(variables = vars, x_train = x_train, y_train = y_train,
                              x_test = x_test, y_test = y_test, arg_lst = arg_lst,
                              missing_threshold = missing_threshold)
            })
            performance_lst[[m]] <- as_tibble(rbindlist(lapply(cvaucs, function(x) {
                tibble::tibble(perf = x$cvAUC, var = x$se ^ 2 * n)
            }))) %>%
                mutate(type = c("gFWER", "PFP", "FDR"))
        } else {
            cvauc <- one_pred_perf(variables = these_selected_vars, x_train = x_train,
                                   y_train = y_train, x_test = x_test,
                                   y_test = y_test, arg_lst = arg_lst,
                                   missing_threshold = missing_threshold)
            performance_lst[[m]] <- tibble::tibble(perf = cvauc$cvAUC, var = cvauc$se ^ 2 * n)
        }
    }
    return(performance_lst)
}

# get prediction performance for a single set of variables
# @param variables the set of variables
# @param x_train the training covariates
# @param x_test the testing covariates
# @param y_train the training outcome
# @param y_test the testing outcome
# @param prediction_procedure the prediction procedure to use (glm, SL)
# @param arg_lst a list of arguments to pass to Super Learner
# @param missing_threshold the missing data threshold
one_pred_perf <- function(variables, x_train, y_train, x_test, y_test,
                          prediction_procedure = "glm", arg_lst,
                          missing_threshold = 0.1) {
    var_nms <- names(x_train)
    selected_variables <- var_nms[variables]
    x_train_selected <- x_train %>%
        select(all_of(selected_variables))
    x_test_selected <- x_test %>%
        select(all_of(selected_variables))
    # train the prediction algorithm
    if (ncol(x_train_selected) == 0) {
        preds <- rep(mean(y_train), length(y_test))
    } else if (prediction_procedure == "glm") {
        fit <- glm(y ~ ., data = data.frame(y = y_train, x_train_selected))
        preds <- predict(fit, newdata = x_test_selected)
    } else {
        if (ncol(x_train_selected) < 2) {
            sl_lib <- c(arg_lst$SL.library[!grepl("glmnet", arg_lst$SL.library)],
                    "SL.glm")
        } else {
            sl_lib <- arg_lst$SL.library
        }
        fit <- SuperLearner::SuperLearner(Y = y_train, X = x_train_selected, SL.library = sl_lib, family = arg_lst$family, method = arg_lst$method, cvControl = arg_lst$cvControl)
        preds <- get_sl_preds(newdata = x_test_selected, fit = fit)
    }
    cvauc <- get_cvauc(preds = preds, y = y_test, missing_threshold = missing_threshold)
    return(cvauc)
}
