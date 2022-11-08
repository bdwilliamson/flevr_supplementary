# get the true performance
# @param test_dat the test dataset
# @param linear whether or not we used a linear conditional mean model
# @param family the type of data ("binomial" implies binary data, otherwise it's continuous)
# @param beta the true value of the coefficient vector
# @param x_dist the x distribution
# @return the true performance
get_true_performance <- function(test_dat, linear = TRUE, family = "binomial", beta, 
                                 x_dist = "normal", y_dist = "many") {
  if (linear) {
    linear_predictor <- cbind(1, as.matrix(test_dat %>% select(-y))) %*% beta
  } else {
    linear_predictor <- nl_conditional_mean(x = test_dat %>% select(-y), beta = beta, 
                                            x_dist = x_dist, y_dist = y_dist)
  }
  if (grepl("binomial", family)) {
    if (grepl("probit", family)) {
      preds <- pnorm(linear_predictor)
    } else {
      preds <- expit(linear_predictor)
    }
    true_perf <- cvAUC::cvAUC(preds, test_dat$y)$cvAUC
  } else {
    true_perf <- vimp::measure_r_squared(linear_predictor, test_dat$y)$point_est
  }
  true_perf
}
