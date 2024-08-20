# Data analysis: investigate cross-sectional VIM performance

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

source(paste0(proj_root, "/code/sims/utils.R"))
source(paste0(proj_root, "/code/data_analysis/00_utils.R"))

# set up args ------------------------------------------------------------------
parser <- OptionParser()
# parser <- add_option(parser, "--run-leave-out", type = "integer", default = 1,
#                      help = "Should we run leave-out groups?")
# parser <- add_option(parser, "--run-add-in", type = "integer", default = 1,
#                      help = "Should we run add-in groups?")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

data_dir <- "G:/CTRHS/IMATS/Data/SRS3 IMATS data/"
results_dir <- "G:/CTRHS/IMATS/Brian/longitudinal_vim/results/data_analysis/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# read in the dataset ----------------------------------------------------------
analysis_dataset <- readRDS(paste0(data_dir, "imats_srs3_cohort_study_analysis_dataset.rds"))
# set up outcome, covariates
y <- analysis_dataset$event90
X <- analysis_dataset %>% 
  select(-study_id, -visit_n, -timepoint, -event90)
timepoints <- unique(analysis_dataset$timepoint)
num_timepoints <- length(timepoints)
# get the people with all 6 visits -- these are the people SEs will be based on
complete_obs <- analysis_dataset %>% 
  group_by(study_id) %>% 
  summarize(num_visits = n()) %>% 
  filter(num_visits == 6) %>% 
  mutate(is_complete = TRUE) %>% 
  select(-num_visits)
is_complete_obs <- analysis_dataset %>% 
  left_join(complete_obs, by = "study_id") %>% 
  mutate(is_complete = as.numeric(ifelse(is.na(is_complete), FALSE, is_complete))) %>% 
  pull(is_complete)
  

# set up the data analysis -----------------------------------------------------
# set up the variable sets
varset_null <- NULL
varset_demographics <- which(
  names(X) %in% c("age", "hispanic", "census_flag", "site_hpi",
    paste0("visit_pc", c(1, 2)),
    paste0("sex_", c("f", "o", "u")), 
    paste0("race_", c("as", "ba", "hp", "in", "mu", "ot", "un")),
    paste0("income_", c("unknown", "lt25", "lt40")),
    paste0("coll_deg_", c("unknown", "lt25")),
    paste0("ins_", c("aca_miss", "aca_num", 
                     "medicaid_miss", "medicaid_num", 
                     "medicare_miss", "medicare_num", 
                     "privatepay_miss", "privatepay_num",
                     "statesubsidized_miss", "statesubsidized_num",
                     "selffunded_miss", "selffunded_num",
                     "highdeductible_miss", "highdeductible_num",
                     "other_miss", "other_num", 
                     "commercial_miss", "commercial_num")),
    "days_since_visit1", "mths_since_prev", "enr_calc")
)
varset_age <- which(names(X) %in% "age")
varset_age_sex <- which(names(X) %in% c("age", paste0("sex_", c("f", "o", "u"))))
varset_age_sex_race_ethnicity <- which(names(X) %in% c("age", paste0("sex_", c("f", "o", "u")),
                                                       paste0("race_", c("as", "ba", "hp", "in", "mu", "ot", "un"))))
varset_prior_self_harm <- which(grepl("asa_", names(X)) | grepl("osa_", names(X)) | 
                                  grepl("lsa_", names(X)) | grepl("aip_", names(X)))
varset_all_diagnoses_utilization <- which(
  (grepl("_visit", names(X)) & !grepl("item9", names(X)) & !grepl("days_since", names(X))) |
    grepl("_days_", names(X)) | grepl("_mths_", names(X))
)
varset_diagnoses <- varset_all_diagnoses_utilization[!(varset_all_diagnoses_utilization %in% varset_prior_self_harm)]
varset_charlson <- which(grepl("charlson_score", names(X)))
varset_phq9 <- which(grepl("item9_visit", names(X)))
varsets <- c(
  list(varset_null), # (1) comparator (for add-in compared to nothing), corresponds to variable set 1 in the paper
  list(varset_demographics), # (2) comparator (for add-in compared to demographics)
  list(c(varset_demographics, varset_prior_self_harm)), # (3) baseline model: demographics + prior self-harm
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses)), # (4) add diagnoses & utilization to baseline model (includes PHQ-8)
  list(c(varset_demographics, varset_prior_self_harm, varset_charlson)), # (5) add charlson to baseline model
  list(c(varset_demographics, varset_prior_self_harm, varset_phq9)), # (6) add PHQ-9 to baseline model
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson)), # (7) add charlson to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_phq9)), # (8) add PHQ-9 to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_charlson, varset_phq9)), # (9) add PHQ-9 to model (5)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson, varset_phq9)), # (10) all variables
  list(varset_phq9), # (11) PHQ-9 compared to nothing, corresponds to variable set 2 in the paper
  list(varset_age), # (12) age compared to nothing
  list(varset_age_sex), # (13) age and sex compared to nothing, corresponds to variable set 3 in the paper
  list(varset_age_sex_race_ethnicity), # (14) age, sex, race and ethnicity compared to nothing
  list(c(varset_age, varset_phq9)), # (15) PHQ-9 compared to age alone (12)
  list(c(varset_age_sex, varset_phq9)), # (16) PHQ-9 compared to age and sex (13), corresponds to variable set 4 in the paper
  list(c(varset_age_sex_race_ethnicity, varset_phq9)), # (17) PHQ-9 compared to age, sex, race and ethnicity (14)
  list(c(varset_age_sex, varset_prior_self_harm)), # (18) corresponds to variable set 5 in the paper
  list(c(varset_age_sex, varset_prior_self_harm, varset_phq9)) # (19) corresponds to variable set 6 in the paper
)
num_varsets <- length(varsets)

# set up the SL library and algorithm
xgb_tune_params <- list(max_depth = c(1, 4), shrinkage = 0.1, ntrees = 500)
rf_tune_params <- list(num.trees = 500, min.node.size = c(1, 10, 50))
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params, 
                               detailed_names = TRUE, name_prefix = "xgb")
rf_learners <- create.Learner("SL.ranger", tune = rf_tune_params, 
                              detailed_names = TRUE, name_prefix = "rf")
learner_lib <- c("SL.mean", "SL.glm", "SL.glmnet", rf_learners$names, xgb_learners$names)
num_learners <- length(learner_lib) + 2
all_learners <- c(learner_lib, "Discrete SL", "SL")
unique_learners <- c(get_unique_learners(learner_lib), "Discrete SL", "SL")
num_unique_learners <- length(unique_learners)

k_outer <- 10
k_inner <- k_outer

sl_opts <- list("family" = "binomial", "method" = "method.CC_nloglik",
                "cvControl" = list(V = k_outer, stratifyCV = TRUE),
                "innerCvControl" = rep(list(list(V = k_inner, stratifyCV = TRUE)), k_outer))
measure_type <- "auc"

# estimate prediction functions, performance for each variable set at each time point -----------------------------
# set up parallelization
num_cores <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(num_cores)
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterExport(cl, learner_lib)

# set up the random number stream
set.seed(20230705)
timepoint_seeds <- round(runif(n = num_timepoints, min = 1e4, max = 1e5))
varset_seeds <- lapply(as.list(1:num_timepoints), function(t) round(runif(n = num_varsets, min = 1e4, max = 1e5)))
  
# set up the return object
output_list <- vector("list", length = num_varsets)

cv_sls <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
cv_folds <- vector("list", length = num_timepoints)
set.seed(20230706)
ss_folds <- lapply(as.list(1:num_timepoints), function(t) vimp::make_folds(1:k_outer, V = 2))
cv_pred_list <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
cv_pred_perf_list <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
cv_vim_list <- vector("list", length = num_timepoints)

# obtain CV SL for each timepoint and varset
start <- Sys.time()
for (i in seq_len(num_timepoints)) {
  this_X <- X[analysis_dataset$timepoint == timepoints[i], ]
  this_y <- y[analysis_dataset$timepoint == timepoints[i]]
  set.seed(timepoint_seeds[i])
  cv_folds[[i]] <- vimp::make_folds(this_y, V = k_outer, stratified = TRUE)
  these_sl_opts <- sl_opts
  for (j in seq_len(num_varsets)) {
    if (!is.null(cv_folds[[i]])) {
      these_sl_opts$cvControl$validRows <- make_cv_sl_folds(cv_folds[[i]])
    }
    this_x_df <- this_X[, varsets[[j]]]
    set.seed(varset_seeds[[i]][j])
    cv_sls[[i]][[j]] <- CV.SuperLearner(
      Y = this_y, X = as.matrix(this_x_df), SL.library = learner_lib,
      family = these_sl_opts$family, method = these_sl_opts$method,
      cvControl = these_sl_opts$cvControl, innerCvControl = these_sl_opts$innerCvControl,
      parallel = cl
    )
    saveRDS(cv_sls[[i]][[j]], file = paste0(results_dir, "cv_sls_", i, "_", j, ".rds"))
  }
  saveRDS(cv_folds, file = paste0(results_dir, "cv_folds.rds"))
}
end <- Sys.time()
cat("Elapsed time: ", format(end - start), "\n")

# compute predictiveness, VIM for each variable set at each time point
for (i in seq_len(num_timepoints)) {
  this_X <- X[analysis_dataset$timepoint == timepoints[i], ]
  this_y <- y[analysis_dataset$timepoint == timepoints[i]]
  these_complete_obs <- is_complete_obs[analysis_dataset$timepoint == timepoints[i]]
  cv_folds_vec <- cv_folds[[i]]
  cv_vim_list[[i]] <- vector("list", length = num_varsets - 1)
  cv_pred_perf_list[[i]][[1]] <- get_all_predictiveness(
    sl_fit = cv_sls[[i]][[1]], type = measure_type,
    complete_obs = these_complete_obs
  )
  # fix for SL.glm not working for null set
  if (nrow(cv_pred_perf_list[[i]][[1]]) != num_learners) {
    tmp <- cv_pred_perf_list[[i]][[1]]
    tmp <- rbind.data.frame(tmp, do.call(rbind.data.frame, rep(list(tmp[tmp$Learner == "SL.mean", ]), num_learners - nrow(tmp))))
    missing_learners <- all_learners[unlist(lapply(as.list(all_learners), function(learner) !any(grepl(learner, cv_pred_perf_list[[i]][[1]]$Learner))))]
    tmp$Learner <- c("SL", "Discrete SL", "SL.mean", missing_learners)  
    cv_pred_perf_list[[i]][[1]] <- tmp
  }
  for (j in 2:num_varsets) {
    this_cv_pred_perf <- get_all_predictiveness(
      sl_fit = cv_sls[[i]][[j]], type = measure_type, complete_obs = these_complete_obs
    )
    cv_pred_perf_list[[i]][[j]] <- this_cv_pred_perf
    full_preds <- get_all_learner_preds(cv_sls[[i]][[j]], learners = unique_learners,
                                        perf_obj = this_cv_pred_perf)
    # if variable set 2: compare to null variable set
    if (j %in% c(2)) {
      comp_indx <- 1
    } else if (j == 3) {
      comp_indx <- 2
    # if variable set 4, 5, 6: compare to baseline
    } else if (j %in% c(4, 5, 6)) {
      comp_indx <- 3
    # if variable set 7, 8: compare to baseline + diagnoses & utilization
    } else if (j %in% c(7, 8)) {
      comp_indx <- 4
    # if variable set 9: compare to baseline + charlson
    } else if (j == 9) {
      comp_indx <- 5
    # if variable set 10: compare to baseline + diagnoses & utilization + charlson
    } else if (j == 10) {
      comp_indx <- 7
    # if variable set 11, compare to null set
    } else if (j %in% c(11, 12)) {
      comp_indx <- 1
    } else if (j %in% c(13, 14)) {
      comp_indx <- j - 1
    } else if (j %in% c(15, 16, 17)) {
      comp_indx <- j - 3
    }
    reduced_predictiveness <- get_all_predictiveness(
      sl_fit = cv_sls[[i]][[comp_indx]], type = measure_type, complete_obs = these_complete_obs
    )
    # fix for SL.glm not working for null set, other learners not working sometimes
    if (nrow(reduced_predictiveness) != num_learners) {
      tmp <- reduced_predictiveness
      tmp <- rbind.data.frame(tmp, do.call(rbind.data.frame, rep(list(tmp[tmp$Learner == "SL.mean", ]), num_learners - nrow(tmp))))
      missing_learners <- all_learners[unlist(lapply(as.list(all_learners), function(learner) !any(grepl(learner, reduced_predictiveness$Learner))))]
      tmp$Learner <- c(reduced_predictiveness$Learner, missing_learners)  
      reduced_predictiveness <- tmp
    }
    reduced_preds <- get_all_learner_preds(cv_sls[[i]][[comp_indx]], 
                                           learners = unique_learners,
                                           perf_obj = reduced_predictiveness)
    all_na <- unlist(lapply(reduced_preds, function(x) all(is.na(x))))
    if (any(all_na)) {
      reduced_preds <- lapply(reduced_preds, function(pred) {
        if (all(is.na(pred))) {
          reduced_preds$mean
        } else {
          pred
        }
      })
    }
    cv_vim_list[[i]][[j - 1]] <- get_all_vims(
      y = this_y, x = this_X, full_preds = full_preds, reduced_preds = reduced_preds,
      var_set = varsets[[j]], measure_type = measure_type, cv_folds = cv_folds_vec,
      ss_folds = ss_folds, alpha = 0.05, K = k_outer / 2, complete_obs = these_complete_obs
    )
  }
}
# reorganize the lists to be in the following order: algorithm, variable set, timepoint
reordered_cv_vim_list <- vector("list", length = num_unique_learners)
for (k in seq_len(num_unique_learners)) {
  reordered_cv_vim_list[[k]] <- vector("list", length = num_varsets - 1)
  for (j in 1:length(cv_vim_list[[1]])) {
    reordered_cv_vim_list[[k]][[j]] <- lapply(cv_vim_list, function(l) l[[j]][[k]])
  }
}
# just need to reorder to variable set, timepoint
reordered_cv_perf_list <- vector("list", length = num_varsets)
for (j in 1:length(cv_pred_perf_list[[1]])) {
  reordered_cv_perf_list[[j]] <- lapply(cv_pred_perf_list, function(l) l[[j]])
}

# estimate longitudinal summaries of predictiveness, VIM for each algo and group of variables
lvim_list <- get_all_lvims(cv_vims = reordered_cv_vim_list, num_timepoints = num_timepoints)

# finalize output and return ---------------------------------------------------
output <- data.table::rbindlist(
  lapply(as.list(seq_len(num_unique_learners)), function(i) {
    collapse_output(output_list = lvim_list[[i]], algo = unique_learners[i],
                        varsets = varsets, vim_types = rep("addi", length(varsets) - 1),
                        baseline_vars = NULL, all_vars = varsets[[10]])
  })
)
saveRDS(output, file = paste0(results_dir, "lvim_output.rds"))
cat("Data analysis complete!\n")
