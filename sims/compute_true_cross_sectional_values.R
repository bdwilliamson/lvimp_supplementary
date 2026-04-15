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
source(here::here("sims", "true_value_utils.R"))
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

metrics <- c("auc", "ppv", "sensitivity")
cutoff_prob <- 0.95


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
                    confounder_beta = confounder_beta, 
                    varset_names = varset_names,
                    cutoff_prob = cutoff_prob, dgm = 1)
  }, future.seed = seeds
)
truths_cor_0_df <- do.call(rbind, truths_cor_0_list)
truths_cor_0 <- truths_cor_0_df %>% 
  group_by(metric, varset, designation, corr_within, corr_between) %>% 
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
                    confounder_beta = confounder_beta, 
                    varset_names = varset_names,
                    cutoff_prob = cutoff_prob, dgm = 1)
  }, future.seed = seeds
)
truths_cor_50_df <- do.call(rbind, truths_cor_50_list)
truths_cor_50 <- truths_cor_50_df %>% 
  group_by(metric, varset, designation, corr_within, corr_between) %>% 
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
                    corr_within = cor_within, confounder_beta = confounder_beta,
                    dgm = 1)
    y_mat <- do.call(cbind, lapply(as.list(1:num_timepoints), function(z) dat$y[dat$t == z]))
    y_cor <- cor(y_mat)
    y_mean <- colMeans(y_mat)
    x_mats <- lapply(as.list(1:p), function(j) {
      do.call(cbind, lapply(as.list(1:num_timepoints), function(z) {
        dat %>% 
          filter(t == z) %>% 
          pull(!!paste0("X", j))
      }))
    })
    x_cor <- lapply(x_mats, cor)
    return(list("y" = y_cor, "x" = x_cor, "y_mean" = y_mean))
  }, future.seed = seeds
)
all_y_cor <- lapply(all_corr_list, function(x) x$y)
all_x_cor <- lapply(as.list(1:10), function(j) {
  lapply(all_corr_list, function(z) z$x[[j]])
})

y_array <- array(as.numeric(unlist(all_y_cor)), dim = c(num_timepoints, num_timepoints, nsim))
y_summ <- apply(y_array, c(1, 2), mean)
saveRDS(y_summ, here::here("..", "results", "sims", "y_cor_dgm_1.rds"))

x_array_list <- lapply(all_x_cor, function(x) array(as.numeric(unlist(x)), dim = c(num_timepoints, num_timepoints, nsim)))
x_summ <- lapply(x_array_list, function(x) apply(x, c(1, 2), mean))
saveRDS(x_summ, here::here("..", "results", "sims", "x_cor_dgm_1.rds"))

y_summ
x_summ

# a new simulation scenario, for studying PPV and sensitivity away from the boundary ---------------
# assoc between confounders and truly important Xs
# const <- 0.075 # AUC, etc. a bit too small / not too different from baseline variables
cor_between <- 0
cor_within <- 0.5

const_dgm2 <- 0.1
beta_01_dgm2 <- rep(const_dgm2, num_timepoints)
beta_02_dgm2 <- const_dgm2 + timepoints / 15
beta_03_dgm2 <- (-1/4) * (1 + exp((-1) * timepoints))^(-1) + const_dgm2
beta_0c_dgm2 <- rep(0.01, num_timepoints)
beta_0_dgm2 <- lapply(as.list(seq_len(num_timepoints)), function(t) {
  matrix(c(beta_01_dgm2[t], beta_02_dgm2[t], beta_03_dgm2[t], rep(beta_0c_dgm2[t], 4), rep(0, p - 7)))
})
lapply(as.list(1:4), function(t) sum(beta_0_dgm2[[t]]))

set.seed(20251002)
x <- gen_x(n = n, p = p, T = num_timepoints, corr_between = cor_between,
           corr_within = cor_within, confounder_beta = confounder_beta, dgm = 2)
lapply(as.list(1:4), function(t) summary(pnorm(as.matrix(x[[t]]) %*% beta_0_dgm2[[t]] - 0.2))) # 0.5 worked well with .075, probs not too small; but with 0.1, too small
cor(x[[1]][,1], x[[2]][,1])
cor(x[[1]][,1], x[[3]][,1])
cor(x[[1]][,1], x[[4]][,1])
# correlations among Xs are the same
set.seed(20250930)
y <- gen_ar_y(n = n, x = x, beta_0 = beta_0_dgm2, T = num_timepoints, dgm = 2, rho = 0.3)
mean(y[[1]])
cor(y[[1]], y[[2]])

set.seed(20250930); tmp <- get_true_values(iteration = i, n = n, p = p, outcome_type = "binary",
                                           num_timepoints = num_timepoints, beta_0 = beta_0_dgm2,
                                           cor_between = cor_between, cor_within = cor_within,
                                           confounder_beta = confounder_beta, 
                                           varset_names = varset_names,
                                           cutoff_prob = cutoff_prob,
                                           dgm = 2)

tmp %>% filter(metric == "auc", varset == "all" | varset == "baseline")
tmp %>% filter(metric == "ppv", varset == "all" | varset == "baseline")
tmp %>% filter(metric == "sensitivity", varset == "all" | varset == "baseline")

cor_between <- 0
cor_within <- 0.5
future::plan(multisession)
current_seed <- 1234
seeds <- future_lapply(as.list(seq_len(nsim)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
truths_cor_50_dgm_2_list <- future.apply::future_lapply(
  X = as.list(seq_len(nsim)), FUN = function(i) {
    get_true_values(iteration = i, n = n, p = p, outcome_type = "binary",
                    num_timepoints = num_timepoints, beta_0 = beta_0_dgm2,
                    cor_between = cor_between, cor_within = cor_within,
                    confounder_beta = confounder_beta, 
                    varset_names = varset_names,
                    cutoff_prob = cutoff_prob, dgm = 2)
  }, future.seed = seeds
)
truths_cor_50_dgm_2_df <- do.call(rbind, truths_cor_50_dgm_2_list)
truths_cor_50_dgm_2 <- truths_cor_50_dgm_2_df %>% 
  group_by(metric, varset, designation, corr_within, corr_between) %>% 
  rename(truth_init = truth) %>% 
  summarize(truth = mean(truth_init), var_truth = var(truth_init), .groups = "drop")
saveRDS(truths_cor_50_dgm_2_df, here::here("..", "results", "sims", "truths_cross_sectional_0_0.5_all_dgm_2.rds"))
saveRDS(truths_cor_50_dgm_2, here::here("..", "results", "sims", "truths_cross_sectional_0_0.5_dgm_2.rds"))

# describe correlation between outcomes, features over time
outcome_corrs <- vector("list", length = nsim)
feature_corrs <- vector("list", length = nsim)
current_seed <- 20240325
future::plan(multisession)
seeds <- future_lapply(as.list(seq_len(nsim)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)
all_corr_list_dgm_2 <- future.apply::future_lapply(
  X = as.list(seq_len(nsim)), FUN = function(i) {
    dat <- gen_data(n = n, p = p, outcome_type = "binary", T = num_timepoints,
                    beta_0 = beta_0_dgm2, corr_between = cor_between, 
                    corr_within = cor_within, confounder_beta = confounder_beta,
                    dgm = 2)
    y_mat <- do.call(cbind, lapply(as.list(1:num_timepoints), function(z) dat$y[dat$t == z]))
    y_cor <- cor(y_mat)
    y_mean <- colMeans(y_mat)
    x_mats <- lapply(as.list(1:p), function(j) {
      do.call(cbind, lapply(as.list(1:num_timepoints), function(z) {
        dat %>% 
          filter(t == z) %>% 
          pull(!!paste0("X", j))
      }))
    })
    x_cor <- lapply(x_mats, cor)
    return(list("y" = y_cor, "x" = x_cor, "y_mean" = y_mean))
  }, future.seed = seeds
)
all_y_means_dgm_2 <- lapply(all_corr_list_dgm_2, function(x) x$y_mean)
all_y_cor_dgm_2 <- lapply(all_corr_list_dgm_2, function(x) x$y)
all_x_cor_dgm_2 <- lapply(as.list(1:10), function(j) {
  lapply(all_corr_list_dgm_2, function(z) z$x[[j]])
})

y_array_dgm_2 <- array(as.numeric(unlist(all_y_cor_dgm_2)), dim = c(num_timepoints, num_timepoints, nsim))
y_summ_dgm_2 <- apply(y_array_dgm_2, c(1, 2), mean)
saveRDS(y_summ_dgm_2, here::here("..", "results", "sims", "y_cor_dgm_2.rds"))

y_means_dgm_2 <- array(as.numeric(unlist(all_y_means_dgm_2)), dim = c(num_timepoints, nsim))
y_means_summ_dgm_2 <- rowMeans(y_means_dgm_2)
saveRDS(y_means_summ_dgm_2, here::here("..", "results", "sims", "y_mean_dgm_2.rds"))


x_array_list_dgm_2 <- lapply(all_x_cor_dgm_2, function(x) array(as.numeric(unlist(x)), dim = c(num_timepoints, num_timepoints, nsim)))
x_summ_dgm_2 <- lapply(x_array_list_dgm_2, function(x) apply(x, c(1, 2), mean))
saveRDS(x_summ_dgm_2, here::here("..", "results", "sims", "x_cor_dgm_2.rds"))

y_summ_dgm_2
y_means_summ_dgm_2
x_summ_dgm_2
