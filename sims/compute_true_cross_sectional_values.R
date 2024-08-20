# Compute true values for the simulation

library("here")
library("dplyr")
library("tidyr")
library("future.apply")

source(here::here("sims", "gen_data.R"))
# set up the data
nsim <- 2500
n <- 1e6
p <- 10
num_timepoints <- 4
timepoints <- seq_len(num_timepoints) - 1
# assoc between confounders and truly important Xs
# confounder_beta <- c(0.5, 0.25, 0.15)
confounder_beta <- c(0.05, 0.05, 0.05, .05)
beta_01 <- rep(2, num_timepoints)
beta_02 <- 2 + timepoints / 4
beta_03 <- (-1) * (1 + exp((-1) * timepoints))^(-1) + 2
beta_0c <- rep(0.05, num_timepoints)
beta_0 <- lapply(as.list(seq_len(num_timepoints)), function(t) {
  matrix(c(beta_01[t], beta_02[t], beta_03[t], rep(beta_0c[t], 4), rep(0, p - 7)))
})

varset_confounders <- 4:7
vars_to_measure <- c(1:3, 8:10)
varsets <- c(list(varset_confounders), 
             lapply(seq_len(length(vars_to_measure)), function(i) {
               sort(c(vars_to_measure[i], varset_confounders))
             }),
             lapply(seq_len(length(vars_to_measure)), function(i) {
               sort((1:10)[-vars_to_measure[i]])
             }))
important_vars <- 1:3
noise_vars <- 8:10
varset_names <- c("all", "baseline", 
                  "addi_1", "addi_2", "addi_3", "addi_8", "addi_9", "addi_10",
                  "loco_1", "loco_2", "loco_3", "loco_8", "loco_9", "loco_10")
varsets_minus_confounders <- varsets[-1]
varsets_minus_confounders_names <- varset_names[-c(1:2)]

# get multiplier for correlation across time
get_multiplier <- function(p, corr, timepoint) {
  if (timepoint == 1) {
    return(diag(p))
  }
  if (timepoint == 2) {
    return(diag(p) + corr * diag(p))
  } else {
    return(diag(p) + corr * get_multiplier(p, corr, timepoint = timepoint - 1))
  }
}

# get the covariance matrix
# @param cor_between correlation between variables at a given timepoint
# @param cor_within correlation within a variable across timepoints
get_sigma <- function(cor_between = 0, cor_within = 0, confounder_beta = c(0.5 - 1e-5, 0.25, 0.15),
                      timepoint = 1) {
  sigma_11 <- matrix(cor_between, nrow = 3, ncol = 3)
  diag(sigma_11) <- 1
  sigma_22 <- matrix(cor_between, nrow = 4, ncol = 4)
  diag(sigma_22) <- 1
  sigma_33 <- diag(nrow = 3, ncol = 3)
  sigma_12 <- t(matrix(confounder_beta, nrow = 4, ncol = 3, byrow = TRUE) %*% sigma_11)
  sigma_21 <- t(sigma_12)
  sigma_31 <- matrix(0, nrow = 3, ncol = 3)
  sigma_13 <- t(sigma_31)
  sigma_32 <- matrix(0, nrow = 3, ncol = 4)
  sigma_23 <- t(sigma_32)
  sigma <- matrix(NA, nrow = 10, ncol = 10)
  sigma[1:3, 1:3] <- sigma_11
  sigma[4:7, 4:7] <- sigma_22
  sigma[8:10, 8:10] <- sigma_33
  sigma[1:3, 4:7] <- sigma_12
  sigma[1:3, 8:10] <- sigma_13
  sigma[4:7, 1:3] <- sigma_21
  sigma[4:7, 8:10] <- sigma_23
  sigma[8:10, 1:3] <- sigma_31
  sigma[8:10, 4:7] <- sigma_32
  sigma_2 <- get_multiplier(p = 10, corr = cor_within, timepoint = timepoint) %*% sigma
  return(sigma)
}

# get part of the covariance matrix for a conditional normal distribution
get_conditional_sigmas <- function(sigma, indices = 4:7) {
  sigma_11 <- sigma[-indices, -indices, drop = FALSE]
  sigma_22 <- sigma[indices, indices, drop = FALSE]
  sigma_12 <- sigma[-indices, indices, drop = FALSE]
  return(list("sigma_11" = sigma_11, "sigma_12" = sigma_12, "sigma_22" = sigma_22))
}

# for AUC
piecewise_linear_estimate <- function(x) {
  if (!is.matrix(x)) {
    indices <- seq_len(length(x))
    x <- matrix(x, nrow = 1)
  }
  indices <- seq_len(ncol(x))
  return(x[, range(indices)[1]] / 2 + x[, range(indices)[2]] / 2 + sum(x[, 2:(range(indices)[2] - 1)]))
}
# for linear trend
U <- cbind(1, matrix(1:num_timepoints))
trend_matrix <- solve(t(U) %*% U) %*% t(U)

# function to get true values using monte-carlo integration
# @param iteration the monte-carlo iteration
# @param n the sample size of the dataset
# @param p the number of features
# @param outcome_type the outcome type (always binary)
# @param num_timepoints the number of timepoints
# @param beta_0 the true regression parameter for all variables
# @param cor_between the correlation between variables at a given timepoint
# @param cor_within the correlation within a variable over time
# @param confounder_beta the regression parameter for confounders
get_true_values <- function(iteration = 1, n = 1000, p = 10, outcome_type = "binary",
                            num_timepoints = 4, beta_0 = rep(list(1:4), num_timepoints),
                            cor_between = 0, cor_within = 0, confounder_beta = rep(1, 4)) {
  dataset_bin <- gen_data(n = n, p = p, outcome_type = outcome_type, T = num_timepoints,
                          beta_0 = beta_0, corr_between = cor_between, 
                          corr_within = cor_within, confounder_beta = confounder_beta)
  # grand mean: all zero (based on combining conditional normal density and marginal normal density)
  grand_mean <- matrix(rep(0, p))
  
  # notation: f_0tj, t in {1,...,4}, j in {1,2,3,{4,5,6,7}}
  true_pred <- tibble::tibble(t = rep(1:4, length(varsets) + 1),
                                  truth = NA,
                                  varset = rep(varset_names, each = num_timepoints),
                                  corr_between = cor_between,
                                  corr_within = cor_within) %>% 
    mutate(designation = paste0("predictiveness-", t)) 
  for (timepoint in 1:num_timepoints) {
    x_t <- dataset_bin %>% 
      filter(t == timepoint) %>% 
      select(-t, -y) %>% 
      as.matrix()
    y_t <- dataset_bin %>% 
      filter(t == timepoint) %>% 
      pull(y)
    sigma <- get_sigma(cor_between = cor_between, cor_within = cor_within, confounder_beta = confounder_beta,
                       timepoint = timepoint)
    f_0t <- pnorm(x_t %*% beta_0[[timepoint]])
    true_pred[true_pred$t == timepoint & true_pred$varset == "all", ]$truth <- cvAUC::AUC(predictions = f_0t, labels = y_t)
    
    # for each set that we're still conditioning on, need
    # a = stuff that we're still conditioning on
    # U = stuff to average out, a random variable
    a_tc <- as.numeric(x_t[, varset_confounders] %*% beta_0[[timepoint]][varset_confounders, , drop = FALSE])
    conf_cond_sigmas <- get_conditional_sigmas(sigma = sigma, indices = varset_confounders)
    conf_sig12_sig22inv <- conf_cond_sigmas$sigma_12 %*% solve(conf_cond_sigmas$sigma_22)
    mu_starc <- as.numeric(t(grand_mean[-varset_confounders] + conf_sig12_sig22inv %*% t(x_t[, varset_confounders])) %*% 
                             beta_0[[timepoint]][-varset_confounders, , drop = FALSE])
    sigma_starc <- as.numeric(t(beta_0[[timepoint]][-varset_confounders, , drop = FALSE]) %*% 
                                (conf_cond_sigmas$sigma_11 - conf_sig12_sig22inv %*% t(conf_cond_sigmas$sigma_12)) %*% 
                                beta_0[[timepoint]][-varset_confounders, , drop = FALSE])
    f_0tc <- pnorm((a_tc + mu_starc) / sqrt(1 + sigma_starc))
    true_pred[true_pred$t == timepoint & true_pred$varset == "baseline", ]$truth <- cvAUC::AUC(predictions = f_0tc, labels = y_t)
    
    for (j in seq_len(length(varsets_minus_confounders))) {
      # this_varset <- c(j, varset_confounders)
      this_varset <- varsets_minus_confounders[[j]]
      a_tj <- as.numeric(x_t[, this_varset] %*% beta_0[[timepoint]][this_varset, , drop = FALSE])
      this_cond_sigmas <- get_conditional_sigmas(sigma = sigma, indices = this_varset)
      this_sig12_sig22inv <- this_cond_sigmas$sigma_12 %*% solve(this_cond_sigmas$sigma_22)
      mu_starj <- as.numeric(t(grand_mean[-this_varset] + this_sig12_sig22inv %*% 
                                 t(x_t[, this_varset])) %*% 
                               beta_0[[timepoint]][-this_varset, , drop = FALSE])
      sigma_starj <- as.numeric(t(beta_0[[timepoint]][-this_varset, , drop = FALSE]) %*% 
                                  (this_cond_sigmas$sigma_11 - this_sig12_sig22inv %*% t(this_cond_sigmas$sigma_12)) %*% 
                                  beta_0[[timepoint]][-this_varset, , drop = FALSE])
      f_0tj <- pnorm((a_tj + mu_starj) / sqrt(1 + sigma_starj))
      true_pred[true_pred$t == timepoint & true_pred$varset == varsets_minus_confounders_names[j], ]$truth <- cvAUC::AUC(predictions = f_0tj, labels = y_t)
      # eval(parse(text = paste0("auc_0", j, "[timepoint] <- cvAUC::AUC(predictions = f_0tj, labels = y_t)")))
    }
  }
  # get longitudinal summaries of prediction performance
  true_lpred_summs <- true_pred %>% 
    group_by(varset, corr_between, corr_within) %>% 
    summarize(`predictiveness-auc` = piecewise_linear_estimate(truth),
              `predictiveness-average` = mean(truth),
              `predictiveness-trend-intercept` = as.numeric(trend_matrix %*% truth)[1],
              `predictiveness-trend-slope` = as.numeric(trend_matrix %*% truth)[2], .groups = "drop")
  true_lpred <- true_pred %>% 
    select(-t) %>% 
    bind_rows(true_lpred_summs %>% 
                pivot_longer(cols = starts_with("predictiveness"), names_to = "designation", values_to = "truth"))
  # get VIMs
  true_vim <- true_pred %>% 
    pivot_wider(names_from = varset, values_from = truth) %>% 
    mutate(across(starts_with("addi_"), ~ .x - baseline, .names = "vim_{.col}"),
           across(starts_with("loco_"), ~ all - .x, .names = "vim_{.col}")) %>% 
    select(-starts_with("addi"), -starts_with("loco"), -all, -baseline) %>% 
    mutate(designation = gsub("predictiveness", "vim", designation)) %>% 
    pivot_longer(cols = starts_with("vim"), names_to = "varset", values_to = "truth") %>% 
    mutate(varset = gsub("vim_", "", varset))
  # get longitudinal summaries of VIMs
  true_lvim_summs <- true_vim %>% 
    group_by(varset, corr_between, corr_within) %>% 
    summarize(`vim-auc` = piecewise_linear_estimate(truth),
              `vim-average` = mean(truth),
              `vim-trend-intercept` = as.numeric(trend_matrix %*% truth)[1],
              `vim-trend-slope` = as.numeric(trend_matrix %*% truth)[2], .groups = "drop")
  true_lvim <- true_vim %>% 
    select(-t) %>% 
    bind_rows(true_lvim_summs %>% 
                pivot_longer(cols = starts_with("vim"), names_to = "designation", values_to = "truth")) %>% 
    mutate(mc_id = iteration, .before = "truth")
  truths <- true_lpred %>%
    mutate(mc_id = iteration, .before = "truth") %>%
    bind_rows(true_lvim)
  return(truths)
}

# No correlation! --------------------------------------------------------------
cor_between <- 0
cor_within <- 0
future::plan(multisession)
current_seed <- 1234
seeds <- future_lapply(as.list(seq_len(nsim)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
truths_cor_0_list <- future.apply::future_lapply(
  X = as.list(seq_len(nsim)), FUN = function(i) {
    get_true_values(iteration = i, n = n, p = p, outcome_type = "binary",
                    num_timepoints = num_timepoints, beta_0 = beta_0,
                    cor_between = cor_between, cor_within = cor_within,
                    confounder_beta = confounder_beta)
  }, future.seed = seeds
)
truths_cor_0_df <- do.call(rbind, truths_cor_0_list)
truths_cor_0 <- truths_cor_0_df %>% 
  group_by(varset, designation, corr_within, corr_between) %>% 
  summarize(truth = mean(truth), var_truth = var(truth), .groups = "drop")
saveRDS(truths_cor_0_df, here::here("..", "results", "sims", "truths_cross_sectional_0_0_all.rds"))
saveRDS(truths_cor_0, here::here("..", "results", "sims", "truths_cross_sectional_0_0.rds"))

# Correlation between features within a time point -----------------------------
# not doing this, since it just changes the value of the VIM

# Correlation within feature across time ---------------------------------------
# nsim <- 500
cor_between <- 0
cor_within <- 0.5
future::plan(multisession)
current_seed <- 1234
seeds <- future_lapply(as.list(seq_len(nsim)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
truths_cor_50_list <- future.apply::future_lapply(
  X = as.list(seq_len(nsim)), FUN = function(i) {
    get_true_values(iteration = i, n = n, p = p, outcome_type = "binary",
                    num_timepoints = num_timepoints, beta_0 = beta_0,
                    cor_between = cor_between, cor_within = cor_within,
                    confounder_beta = confounder_beta)
  }, future.seed = seeds
)
truths_cor_50_df <- do.call(rbind, truths_cor_50_list)
truths_cor_50 <- truths_cor_50_df %>% 
  group_by(varset, designation, corr_within, corr_between) %>% 
  summarize(truth = mean(truth), var_truth = var(truth), .groups = "drop")
saveRDS(truths_cor_50_df, here::here("..", "results", "sims", "truths_cross_sectional_0_0.5_all.rds"))
saveRDS(truths_cor_50, here::here("..", "results", "sims", "truths_cross_sectional_0_0.5.rds"))

# describe correlation between outcomes, features over time
outcome_corrs <- vector("list", length = nsim)
feature_corrs <- vector("list", length = nsim)
current_seed <- 20240325
future::plan(multisession)
seeds <- future_lapply(as.list(seq_len(nsim)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
all_corr_list <- future.apply::future_lapply(
  X = as.list(seq_len(nsim)), FUN = function(i) {
    dat <- gen_data(n = n, p = p, outcome_type = "binary", T = num_timepoints,
                    beta_0 = beta_0, corr_between = cor_between, 
                    corr_within = cor_within, confounder_beta = confounder_beta)
    y_mat <- do.call(cbind, lapply(as.list(1:num_timepoints), function(z) dat$y[dat$t == z]))
    y_cor <- cor(y_mat)
    x_mats <- lapply(as.list(1:p), function(j) {
      do.call(cbind, lapply(as.list(1:num_timepoints), function(z) {
        dat %>% 
          filter(t == z) %>% 
          pull(!!paste0("X", j))
      }))
    })
    x_cor <- lapply(x_mats, cor)
    return(list("y" = y_cor, "x" = x_cor))
  }, future.seed = seeds
)
all_y_cor <- lapply(all_corr_list, function(x) x$y)
all_x_cor <- lapply(as.list(1:10), function(j) {
  lapply(all_corr_list, function(z) z$x[[j]])
})

y_array <- array(as.numeric(unlist(all_y_cor)), dim = c(num_timepoints, num_timepoints, nsim))
y_summ <- apply(y_array, c(1, 2), mean)

x_array_list <- lapply(all_x_cor, function(x) array(as.numeric(unlist(x)), dim = c(num_timepoints, num_timepoints, nsim)))
x_summ <- lapply(x_array_list, function(x) apply(x, c(1, 2), mean))

y_summ
x_summ
