# Simulation: investigate cross-sectional VIM performance

# load required packages and functions -----------------------------------------
library("tidyverse")
library("data.table")
library("SuperLearner")
library("glmnet")
library("ranger")
library("xgboost")
library("vimp")
library("lvimp") # if you need to update: pak::pkg_install("bdwilliamson/lvimp")
library("parallel")
library("optparse")
library("here")
library("rprojroot")
this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

# source(here::here("sims", "gen_data.R"))
source(paste0(proj_root, "/code/sims/gen_data.R"))
source(paste0(proj_root, "/code/sims/investigate_cross_sectional_performance_once.R"))
source(paste0(proj_root, "/code/sims/utils.R"))
# source(here::here("sims", "investigate_cross_sectional_performance_once.R"))
# source(here::here("sims", "utils.R"))
# set up args ------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--outcome-type", default = "binary",
                     help = "The outcome type (binary or continuous)")
parser <- add_option(parser, "--cor-between", type = "numeric", default = 0,
                     help = "Between-feature correlation")
parser <- add_option(parser, "--cor-within", 
                     type = "numeric", default = 0, 
                     help = "Within-feature correlation")
parser <- add_option(parser, "--n", type = "integer", 
                     default = 500, help = "The sample size")
parser <- add_option(parser, "--p", type = "integer", 
                     default = 10, help = "The number of features")
parser <- add_option(parser, "--num-timepoints", type = "integer", default = 4,
                     help = "Number of timepoints")
parser <- add_option(parser, "--nreps-total", type = "integer", default = 1000,
                     help = "The total number of replications")
parser <- add_option(parser, "--nreps-per-job", type = "integer", default = 1000,
                     help = "The number of replications per job")
parser <- add_option(parser, "--parallel-cv", type = "integer", default = 1,
                     help = "Parallelize the CV.SL (1) or ranger/xgboost (0)?")
parser <- add_option(parser, "--run-leave-out", type = "integer", default = 1,
                     help = "Should we run leave-out groups?")
parser <- add_option(parser, "--run-add-in", type = "integer", default = 1,
                     help = "Should we run add-in groups?")
parser <- add_option(parser, "--simple-model", default = 0, 
                     help = "Should we run simple procedures only?")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

# output_dir <- here::here("..", "..", "results", "sims")
output_dir <- paste0(proj_root, "/results/sims/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# set up the simulation --------------------------------------------------------
all_analyses <- expand.grid(# n = c(5e2, 1e3, 5e3, 1e4),
  n = c(1e2, 2.5e2, 5e2, 1e3, 5e3),
  # p = c(10, 50, 100, 200),
  p = 10,
  cor_between = c(0, 0.2, 0.5),
  cor_within = c(0, 0.2, 0.5),
  outcome_type = c("binary", "continuous"))
n_jobs <- args$nreps_total / args$nreps_per_job
job_id <- which((args$n == all_analyses$n) & (args$p == all_analyses$p) &
                  (args$cor_between == all_analyses$cor_between) &
                  (args$cor_within == all_analyses$cor_within) &
                  (args$outcome_type == all_analyses$outcome_type))
save_iter <- 250
if (args$n == 1e4) {
  job_id <- which((args$p == all_analyses$p) &
                    (args$cor_between == all_analyses$cor_between) &
                    (args$cor_within == all_analyses$cor_within) &
                    (args$outcome_type == all_analyses$outcome_type))[1] + 1e4
  save_iter <- 50
}
# set up the effect sizes
timepoints <- seq_len(args$num_timepoints) - 1
confounder_beta <- c(0.05, 0.05, 0.05, .05)
beta_01 <- rep(2, args$num_timepoints)
beta_02 <- 2 + timepoints / 4
beta_03 <- (-1) * (1 + exp((-1) * timepoints))^(-1) + 2
beta_0c <- rep(0.05, args$num_timepoints)
beta_0 <- lapply(as.list(seq_len(args$num_timepoints)), function(t) {
  matrix(c(beta_01[t], beta_02[t], beta_03[t], rep(beta_0c[t], 4), rep(0, args$p - 7)))
})

# set up the variable sets
varset_confounders <- 4:7
varsets <- NULL
vim_types <- NULL
vars_to_measure <- c(1:3, 8:10)
# add-in variable importance
if (args$run_add_in) {
  varsets <- c(list(varset_confounders),
               lapply(seq_len(length(vars_to_measure)), function(i) {
                 sort(c(vars_to_measure[i], varset_confounders))
               }))
  vim_types <- rep(list("addi"), length(varsets))
}
if (args$run_leave_out) {
  loco_varsets <- c(list(1:10),
                    lapply(seq_len(length(vars_to_measure)), function(i) {
                      sort((1:10)[-vars_to_measure[i]])
                    }))
  varsets <- c(varsets, loco_varsets)
  vim_types <- c(vim_types, rep(list("loco"), length(loco_varsets)))
}

# set up the SL library
num_cores <- parallel::detectCores()
xgb_tune_params <- list(max_depth = 1, shrinkage = 0.1, ntrees = 500)
# xgb_tune_params <- list(max_depth = 4, shrinkage = c(1e-3, 1e-2, 1e-1), ntrees = 500)
rf_tune_params <- list(num.trees = 500)
# rf_tune_params <- list(num.trees = 500, min.node.size = c(1, 10, 25, 50))
if (!as.logical(args$parallel_cv)) {
  xgb_tune_params <- c(xgb_tune_params, list(nthread = num_cores))
  rf_tune_params <- c(rf_tune_params, list(num.threads = num_cores))
}
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params, 
                               detailed_names = TRUE, name_prefix = "xgb")
rf_learners <- create.Learner("SL.ranger", tune = rf_tune_params, 
                              detailed_names = TRUE, name_prefix = "rf")
learner_lib <- c("SL.glm", "SL.glmnet", rf_learners$names, xgb_learners$names)
if (args$n >= 1e4) {
  learner_lib <- c("SL.glm", "SL.glmnet")
}
if (args$simple_model == 1) {
  learner_lib <- "SL.glm"
}

# run the simulation the specified number of times -----------------------------
k_outer <- ifelse(args$n > 1000, 5, 10)
k_inner <- k_outer
# k_inner <- ifelse(args$n <= 500, 20, 10)
if (args$n >= 1e4) {
  k_outer <- 5
  k_inner <- 5
}
# set up parallelization
cl <- parallel::makePSOCKcluster(num_cores)
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterExport(cl, learner_lib)
# set up the random number stream
seed <- job_id * 1000 # need to add something for intermediate jobs
set.seed(seed)
clusterSetRNGStream(cl = cl, iseed = seed)
# run the simulation
cat("Running analysis\n")
start <- Sys.time()
output_list <- vector("list", length = args$nreps_per_job)
for (i in seq_len(args$nreps_per_job)) {
  output_list[[i]] <- investigate_cross_sectional_performance_once(
    mc_id = i + (n_jobs - 1) * args$nreps_per_job, n = args$n, p = args$p, outcome_type = args$outcome_type, T = args$num_timepoints,
    beta_0 = beta_0, confounder_beta = confounder_beta, 
    corr_between = args$cor_between, corr_within = args$cor_within,
    k_outer = k_outer, k_inner = k_inner,
    learners = learner_lib, varsets = varsets, vim_types = vim_types,
    parallel = as.logical(args$parallel_cv), cl = cl
  )
  if ((i %% save_iter) == 0) {
    cat("Finished iteration ", i, "\n")
    saveRDS(output_list, file = paste0(output_dir, "/output_", args$outcome_type,
                                       "_n_", args$n, "_p_", args$p,
                                       "_cb_", args$cor_between,
                                       "_cw_", args$cor_within,
                                       "_ai_", args$run_add_in,
                                       "_lo_", args$run_leave_out, "_interim.rds"))
  }
}
end <- Sys.time()
cat("Elapsed time: ", format(end - start), "\n")
output <- data.table::rbindlist(output_list)
saveRDS(output, file = paste0(output_dir, "/output_", args$outcome_type,
                              "_n_", args$n, "_p_", args$p,
                              "_cb_", args$cor_between,
                              "_cw_", args$cor_within, 
                              "_ai_", args$run_add_in,
                              "_lo_", args$run_leave_out,
                              "_id_", job_id, ".rds"))
cat("Simulation complete!\n")