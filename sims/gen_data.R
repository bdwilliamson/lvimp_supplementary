# generate data for the cross-sectional variable importance simulations

#' @title Generate data for cross-sectional variable importance simulations
#' 
#' @param n the sample size
#' @param p the number of features
#' @param outcome_type the type of outcome ("continuous" or "binary")
#' @param T the number of timepoints
#' @param beta_0 a list of matrices; the coefficients at each time point
#' @param corr_between the between-variable correlation at a given time point
#' @param corr_within the within-variable correlation across time points
#' @param outcome_corr_type type of outcome correlation (defaults to "none", can be "ar" or "glmm")
#' @param confounder_beta the association between the confounders and the important variables
#' @return a tibble with the dataset
gen_data <- function(n = 100, p = 10, outcome_type = "binary", T = 4,
                     beta_0 = lapply(as.list(seq_len(T)), function(t) {
                      matrix(rep(1, p))
                     }), 
                     corr_between = 0, corr_within = 0, 
                     outcome_corr_type = "none",
                     confounder_beta = c(0.5, 0.25, 0.15)) {
  # generate X -- this function returns a list of x matrices
  x <- gen_x(n = n, p = p, T = T, corr_between = corr_between, corr_within = corr_within,
             confounder_beta = confounder_beta)
  # generate Y
  if (grepl("glmm", outcome_corr_type, ignore.case = TRUE)) {
    # Y is based on a (G)LMM
    if (grepl("strong", outcome_corr_type)) {
      cor <- 0.9
    } else {
      cor <- 0.2
    }
    corMat <- matrix(cor, nrow = T, ncol = T)
    diag(corMat) <- 1
    b <- rnorm(n, mean = 0, sd = 1)
    mu <- lapply(seq_len(T), function(t) {
      b + as.matrix(x[[t]]) %*% beta_0[[t]]
    })
    mu_mat <- do.call(cbind, mu)
    latent_y <- mu_mat + mvtnorm::rmvnorm(n, rep(0, T), corMat)
    y_mat <- apply(latent_y, 2, function(z) as.numeric(z > 0))
    y <- lapply(seq_len(T), function(t) y_mat[, t])
  } else if (grepl("ar", outcome_corr_type, ignore.case = TRUE)) {
    # Y is based on an AR(1) model
    y <- vector("list", length = T)
    y[[1]] <- gen_timepoint_y(x = x[[1]], beta_0 = beta_0[[1]], outcome_type = outcome_type)
    for (t in 2:T) {
      y[[t]] <- gen_timepoint_y(x = x[[t]], beta_0 = beta_0[[t]], outcome_type = outcome_type,
                                outcome_corr_type = outcome_corr_type, previous_x = x[[t - 1]])
    }
  } else {
    # Y is independent across time
    y <- lapply(seq_len(T), function(t) {
      gen_timepoint_y(x = x[[t]], beta_0 = beta_0[[t]], outcome_type = outcome_type)
    })  
  }
  
  # return a "long" dataset
  dataset <- tibble::as_tibble(data.table::rbindlist(
    lapply(seq_len(T), function(t) {
      cbind("t" = t, "y" = y[[t]], x[[t]])
    })
  ))
  return(dataset)
}

#' @title Generate outcomes for a specific time point, given covariates
#' 
#' @param x the time-specific features
#' @param outcome_type the type of outcome ("continuous" or "binary")
#' @param outcome_corr_type type of outcome correlation (defaults to "none", can be "ar" or "glmm")
#' @param b a vector of random intercepts for GLMM correlation
#' @param previous_x covariates from the previous timepoint for AR(1) correlation
#' @return the time-specific outcomes
gen_timepoint_y <- function(x = replicate(10, rnorm(100, 0, 1)), 
                            beta_0 = matrix(c(1, 1, 0.5, rep(.25, 4), 
                                              rep(0, ncol(x) - 7))), 
                            outcome_type = "binary",
                            outcome_corr_type = "none", b = rep(0, nrow(x)),
                            previous_x = replicate(10, rnorm(100, 0, 1))) {
  if (grepl("glmm", outcome_corr_type, ignore.case = TRUE)) {
    linear_predictor <- b + as.matrix(x) %*% beta_0
  } else if (grepl("ar", outcome_corr_type, ignore.case = TRUE)) {
    linear_predictor <- as.matrix(previous_x) %*% beta_0
  } else {
    linear_predictor <- as.matrix(x) %*% beta_0
  }
  noisy_predictor <- linear_predictor + rnorm(n = nrow(x), mean = 0, sd = 1)  
  if (outcome_type == "binary") {
    y <- as.numeric(noisy_predictor > 0)
  } else {
    y <- noisy_predictor
  }
  return(y)
}

#' @title Generate errors for covariates
#' 
#' @param n the number of observations
#' @param mean the mean values to generate from
#' @param r_between the between-variable correlation (cholesky decomposition of a covariance matrix)
#' @param p the number of variables to generate
gen_errors <- function(n = 100, mean = rep(0, n), r_between = chol(matrix(0, nrow = p, ncol = p)), p = 3) {
  return(data.frame(matrix(rnorm(n * p, mean = mean, sd = 1), nrow = n, ncol = p) %*% r_between))
}

#' @title Generate features for a specified number of timepoints
#' 
#' @param n the sample size
#' @param p the number of features
#' @param T the number of timepoints
#' @param corr_between the between-variable correlation at a given timepoint
#' @param corr_within the within-variable correlation across timepoints
#' @param confounder_beta the association between the confounders and the important variables
#' @return the time-specific features
gen_x <- function(n = 100, p = 10, T = 4, corr_between = 0, corr_within = 0,
                  confounder_beta = rep(0.05, 4)) {
  # generate T matrices filled with noise variables (note variables 8:p are independent, no importance)
  x <- lapply(as.list(seq_len(T)), function(t) {
    data.frame(matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p))
  })
  # generate (X_4, X_5, X_6, X_7)_1 (confounders)
  sigma_between_confounders <- matrix(corr_between ^ 2, nrow = 4, ncol = 4)
  diag(sigma_between_confounders) <- 1
  r_between_confounders <- chol(sigma_between_confounders)
  x[[1]][, 4:7] <- gen_errors(n = n, mean = rep(0, n), r_between = r_between_confounders, p = 4)
  # generate (X_1, X_2, X_3)_1 (truly important variables)
  imp_means_1 <- as.matrix(x[[1]][, 4:7]) %*% as.matrix(confounder_beta)
  sigma_between_important <- matrix(corr_between ^ 2, nrow = 3, ncol = 3)
  diag(sigma_between_important) <- 1
  r_between_important <- chol(sigma_between_important)
  x[[1]][, 1:3] <- gen_errors(n = n, mean = imp_means_1, r_between = r_between_important, p = 3)
  # generate X for the remaining timepoints
  for (t in seq_len(T)[-1]) {
    x[[t]][, 4:7] <- corr_within * x[[t - 1]][, 4:7] + 
                      gen_errors(n = n, mean = rep(0, n), r_between = r_between_confounders, p = 4)
    imp_means_t <- as.matrix(x[[t]][, 4:7]) %*% as.matrix(confounder_beta)
    x[[t]][, 1:3] <- corr_within * x[[t - 1]][, 1:3] + 
                      gen_errors(n = n, mean = imp_means_t, r_between = r_between_important, p = 3)
  }
  return(x)
}