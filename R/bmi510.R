#' Sample Wrapper for Vectors and DataFrames
#' 
#' A wrapper around sample that tests whether x is an atomic vector or dataframe-like object 
#' and then returns either n samples or n rows as appropriate
#'
#' @param x The object you want to sample from
#' @param n The number of samples you want to draw
#' @param replace Whether to draw samples with or without replacement (default is true)
#' @return The random samples drawn from x
#' @export
#' @examples
#' test_vec = c(5, 6, 7, 8)
#' rando(test_vec, 4, FALSE)
rando = function(x, n = 1, replace = TRUE) {
  calc_sample = NULL
  if(is.data.frame(x)) {
    random_selection = sample(x = nrow(x), size = n, replace = replace)
    calc_sample = x[random_selection,]
    rownames(calc_sample) <- NULL
  } else {
    calc_sample = sample(x = x, size = n, replace = replace)
  }
  return(calc_sample)
}

#' Minimum Value Finder
#' 
#' Accepts an atomic vector x and returns a logical with TRUE where x equals its minimum value.
#'
#' @param x An atomic vector of numbers
#' @param na.rm Whether to omit NAs or not (default is true)
#' @return A logical vector with TRUE being the minimum value
#' @export
#' @examples
#' test_vec = c(7, 6, 8, 5)
#' is_min(test_vec)
is_min = function(x, na.rm = T) {
  # According to which.min/max doc, missing and NaN values are discarded, so this step is redundant
  if(na.rm == T) {
    x = na.omit(x)
  }
  min_index = which.min(x = x)
  final_vec = rep(FALSE, length(x))
  final_vec[min_index] = TRUE
  return(final_vec)
}

#' Maximum Value Finder
#' 
#' Accepts an atomic vector x and returns a logical with 
#' TRUE where x equals its maximum value.
#'
#' @param x An atomic vector of numbers
#' @param na.rm Whether to omit NAs or not (default is true)
#' @return A logical vector with TRUE being the maximum value
#' @export
#' @examples
#' test_vec = c(7, 5, 8, 6)
#' is_max(test_vec)
is_max = function(x, na.rm = T) {
  # According to which.min/max doc, missing and NaN values are discarded, so this step is redundant
  if(na.rm == T) {
    x = na.omit(x)
  }
  max_index = which.max(x = x)
  final_vec = rep(FALSE, length(x))
  final_vec[max_index] = TRUE
  return(final_vec)
}

#' Matrix Repetition from MATLAB
#' 
#' Accepts a dataframe or matrix and x and returns a matrix created by 
#' replicating the rows (or columns) M (N) times.
#'
#' @param x A dataframe or matrix to replicate
#' @param M How many times to replicate the dataframe or matrix row-wise
#' @param N How many times to replicate the dataframe or matrix column-wise
#' @export
#' @examples
#' test_df = data.frame(x = c(1, 2, 3, 4), y = c(5, 6, 7, 8))
#' rep_mat(test_df, 2, 4)
rep_mat = function(x, M = 1, N = 1) {
  new_matrix = x
  if(M != 1) {
    for(row_copy in 1:(M - 1)) {
      new_matrix = rbind(new_matrix, x)
    }
  }
  x = new_matrix
  if(N != 1) {
    for(col_copy in 1:(N - 1)) {
      new_matrix = cbind(new_matrix, x)
    }
  }
  return(new_matrix)
}

#' Classes of Tibble Columns
#' 
#' Returns a character vector containing the classes of each variable in the tibble x.
#'
#' @param x A tibble
#' @export
#' @examples
#' test_tibble = as_tibble(data.frame(a = 1:3, b = letters[1:3], c = Sys.Date() - 1:3))
#' classes(test_tibble)
classes = function(x) {
  base_vector = lapply(x, class)
  fixed_vector = unname(unlist(base_vector))
  return(fixed_vector)
}

#' Numeric Column Scaler
#' 
#' Scales the numeric columns in a tibble.
#'
#' @param x A tibble
#' @export
#' @examples
#' test_tibble = as_tibble(data.frame(a = 1:3, b = letters[1:3], c = Sys.Date() - 1:3))
#' classes(test_tibble)
df_scale = function(x, center = T, scale = T) {
  col_classes = classes(x)
  for(col_index in 1:ncol(x)) {
    if(col_classes[col_index] == "integer" | col_classes[col_index] == "numeric") {
      x[,col_index] = scale(x[,col_index], center = center, scale = scale)
    }
  }
  return(x)
}

#' Normal Log-Likelihood
#' 
#' Calculates the log likelihood of a given dataset belonging to a normal distribution
#'
#' @param x A dataset you want to estimate the likelihood of
#' @param mean The mean of the target normal distribution
#' @param sd The standard deviation of the target normal distribution
#' @export
#' @examples
#' test_norm = rnorm(1000, 4, 2.5)
#' log_likelihood_norm(test_norm, 4, 2.5)
log_likelihood_norm = function(x, mean, sd) {
  result = sum(log(dnorm(x = x, mean = mean, sd = sd)))
  return(result)
}

#' Uniform Log-Likelihood
#' 
#' Calculates the log likelihood of a given dataset belonging to a uniform distribution
#'
#' @param x A dataset you want to estimate the likelihood of
#' @param min The minimum of the target uniform distribution
#' @param max The maximum of the target uniform distribution
#' @export
#' @examples
#' test_unif = runif(1000, 0, 1)
#' log_likelihood_unif(test_unif, 0, 1)
log_likelihood_unif = function(x, min, max) {
  result = sum(log(dunif(x = x, min = min, max = max)))
  return(result)
}

#' Chi-Square Log-Likelihood
#' 
#' Calculates the log likelihood of a given dataset belonging to a chi-squared distribution
#'
#' @param x A dataset you want to estimate the likelihood of
#' @param df The degrees of freedom of the target chi-square distribution
#' @export
#' @examples
#' test_chisq = rchisq(1000, 3)
#' log_likelihood_chisq(test_chisq, 3)
log_likelihood_chisq = function(x, df) {
  result = sum(log(dchisq(x = x, df = df)))
  return(result)
}

#' F Log-Likelihood
#' 
#' Calculates the log likelihood of a given dataset belonging to an F distribution
#'
#' @param x A dataset you want to estimate the likelihood of
#' @param df1 The numerator degrees of freedom of the target F distribution
#' @param df2 The denominator degrees of freedom of the target F distribution
#' @export
#' @examples
#' test_f = rf(1000, 12, 17)
#' log_likelihood_f(test_f, 12, 17)
log_likelihood_f = function(x, df1, df2) {
  result = sum(log(df(x = x, df1 = df1, df2 = df2)))
  return(result)
}

#' T Log-Likelihood
#' 
#' Calculates the log likelihood of a given dataset belonging to a t distribution
#'
#' @param x A dataset you want to estimate the likelihood of
#' @param df The degrees of freedom of the target t distribution
#' @export
#' @examples
#' test_t = rt(1000, 19)
#' log_likelihood_t(test_t, 19)
log_likelihood_t = function(x, df) {
  result = sum(log(dt(x = x, df = df)))
  return(result)
}

#' Sensitivity of a Model
#' 
#' Calculates the sensitivity of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' sensitivity(predicted_vals, ground_truth)
sensitivity = function(pred, truth) {
  true_positives = sum(pred == truth & pred == 1)
  false_negatives = sum(pred != truth & pred == 0)
  calc_sens = true_positives / (true_positives + false_negatives)
  return(calc_sens)
}

#' Specificity of a Model
#' 
#' Calculates the specificity of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' specificity(predicted_vals, ground_truth)
specificity = function(pred, truth) {
  true_negatives = sum(pred == truth & pred == 0)
  false_positives = sum(pred != truth & pred == 1)
  calc_spec = true_negatives / (true_negatives + false_positives)
  return(calc_spec)
}

#' Precision of a Model
#' 
#' Calculates the precision of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' precision(predicted_vals, ground_truth)
precision = function(pred, truth) {
  true_positives = sum(pred == truth & pred == 1)
  false_positives = sum(pred != truth & pred == 1)
  calc_prec = true_positives / (true_positives + false_positives)
  return(calc_prec)
}

#' Recall of a Model
#' 
#' Calculates the recall of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' recall(predicted_vals, ground_truth)
recall = function(pred, truth) {
  true_positives = sum(pred == truth & pred == 1)
  false_negatives = sum(pred != truth & pred == 0)
  calc_rec = true_positives / (true_positives + false_negatives)
  return(calc_rec)
}

#' Accuracy of a Model
#' 
#' Calculates the accuracy of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' accuracy(predicted_vals, ground_truth)
accuracy = function(pred, truth) {
  calc_acc = sum(pred == truth) / length(pred)
  return(calc_acc)
}

#' F1 Score of a Model
#' 
#' Calculates the F1 score of a binary predictor using ground truth and prediction vectors
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
#' predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
#' f1(predicted_vals, ground_truth)
f1 = function(pred, truth) {
  calc_f1 = 2 * precision(pred, truth) * recall(pred, truth) / (precision(pred, truth) + recall(pred, truth))
  return(calc_f1)
}

#' Minimum Sample Size from Cohen's d
#' 
#' Calculates the minimum sample size required to perform a two-sample 
#' t-test with a desired Cohen's effect size.
#'
#' @param d The target effect size
#' @param power The desired power of the statistical test
#' @export
#' @examples
#' minimum_n_per_group(0.76)
minimum_n_per_group = function(d, power = 0.8) {
  power_calculation = power.t.test(delta = d, power = power)
  return(round(power_calculation$n))
}

#' R-Squared Statistic
#' 
#' Calculates the R^2 value for a set of estimated and true quantitative values
#'
#' @param pred The estimated values from the model
#' @param truth The ground truth values
#' @export
#' @examples
#' pred_r2 = c(5.5, 7.7, 9.9, 11.11, 13.13, 15.15)
#' truth_r2 = c(5, 8, 12, 4, 18, 20)
#' r2(pred_r2, truth_r2)
r2 = function(pred, truth) {
  r2_cor = cor(pred, truth) ^ 2
  return(r2_cor)
}

#' Adjusted R-Squared Statistic
#' 
#' Calculates the adjusted R^2 value for a set of estimated and true quantitative values
#'
#' @param pred The predicted values from the model
#' @param truth The ground truth values
#' @param n_p The number of parameters in the model (excluding the intercept)
#' @export
#' @examples
#' pred_r2 = c(5.5, 7.7, 9.9, 11.11, 13.13, 15.15)
#' truth_r2 = c(5, 8, 12, 4, 18, 20)
#' adj_R2(pred_r2, truth_r2, 2)
adj_R2 = function(pred, truth, n_p) {
  sample_size = length(pred)
  adj_r_squared = 1 - ((1 - r2(pred, truth)) * (sample_size - 1)) / (sample_size - n_p - 1)
  return(adj_r_squared)
}