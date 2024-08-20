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
library("doParallel")
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
  list(varset_null), # (1) comparator (for add-in compared to nothing)
  list(varset_demographics), # (2) comparator (for add-in compared to demographics)
  list(c(varset_demographics, varset_prior_self_harm)), # (3) baseline model: demographics + prior self-harm
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses)), # (4) add diagnoses & utilization to baseline model (includes PHQ-8)
  list(c(varset_demographics, varset_prior_self_harm, varset_charlson)), # (5) add charlson to baseline model
  list(c(varset_demographics, varset_prior_self_harm, varset_phq9)), # (6) add PHQ-9 to baseline model
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson)), # (7) add charlson to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_phq9)), # (8) add PHQ-9 to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_charlson, varset_phq9)), # (9) add PHQ-9 to model (5)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson, varset_phq9)), # (10) all variables
  list(varset_phq9), # (11) PHQ-9 compared to nothing
  list(varset_age), # (12) age compared to nothing
  list(varset_age_sex), # (13) age and sex compared to nothing
  list(varset_age_sex_race_ethnicity), # (14) age, sex, race and ethnicity compared to nothing
  list(c(varset_age, varset_phq9)), # (15) PHQ-9 compared to age alone (12)
  list(c(varset_age_sex, varset_phq9)), # (16) PHQ-9 compared to age and sex (13)
  list(c(varset_age_sex_race_ethnicity, varset_phq9)) # (17) PHQ-9 compared to age, sex, race and ethnicity (14)
)
num_varsets <- length(varsets)

# set up the SL library and algorithm
xgb_tune_params <- list(max_depth = c(1, 4), shrinkage = 0.1, ntrees = 500)
rf_tune_params <- list(num.trees = 500, min.node.size = c(1, 10, 50))
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params, 
                               detailed_names = TRUE, name_prefix = "xgb")
rf_learners <- create.Learner("SL.ranger", tune = rf_tune_params, 
                              detailed_names = TRUE, name_prefix = "rf")
# learner_lib <- c("SL.mean", "SL.glm", "SL.glmnet", rf_learners$names, xgb_learners$names)
learner_lib <- c("SL.mean", "SL.glm", "SL.glmnet")
num_learners <- length(learner_lib) + 2
all_learners <- c(learner_lib, "Discrete SL", "SL")
unique_learners <- c(get_unique_learners(learner_lib), "Discrete SL", "SL")
num_unique_learners <- length(unique_learners)

k_outer <- 10
k_inner <- k_outer

sl_opts <- list("family" = "binomial", "method" = "method.CC_nloglik",
                "control" = list(saveFitLibrary = TRUE),
                "cvControl" = list(V = k_outer, stratifyCV = TRUE),
                "innerCvControl" = rep(list(list(V = k_inner, stratifyCV = TRUE)), k_outer))
measure_type <- "auc"

# estimate prediction functions, performance for each variable set at each time point -----------------------------
# set up parallelization
num_cores <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(num_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterExport(cl, learner_lib)

# set up the random number stream
set.seed(20230705)
timepoint_seeds <- round(runif(n = num_timepoints, min = 1e4, max = 1e5))
varset_seeds <- lapply(as.list(1:num_timepoints), function(t) round(runif(n = num_varsets, min = 1e4, max = 1e5)))

# set up the return object
output_list <- vector("list", length = num_varsets)

glms <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
glm_coefs <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
lassos <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
lasso_coefs <- lapply(as.list(1:num_timepoints), function(t) vector("list", num_varsets))
set.seed(20230706)

# obtain lasso, glm for each timepoint and varset
start <- Sys.time()
for (i in seq_len(num_timepoints)) {
  this_X <- X[analysis_dataset$timepoint == timepoints[i], ]
  this_y <- y[analysis_dataset$timepoint == timepoints[i]]
  set.seed(timepoint_seeds[i])
  for (j in seq_len(num_varsets)) {
    this_x_df <- this_X[, varsets[[j]]]
    this_df <- data.frame(y = this_y, this_x_df)
    set.seed(varset_seeds[[i]][j])
    glms[[i]][[j]] <- glm(y ~ ., family = "binomial", data = this_df)
    these_glm_coefs <- coefficients(summary(glms[[i]][[j]]))
    glm_coefs[[i]][[j]] <- as_tibble(these_glm_coefs) |> 
      mutate(Timepoint = i, Varset = j, Variable = rownames(these_glm_coefs))
    if (length(varsets[[j]]) < 3) {
      lassos[[i]][[j]] <- glms[[i]][[j]]
      lasso_coefs[[i]][[j]] <- glm_coefs[[i]][[j]]
    } else {
      lassos[[i]][[j]] <- glmnet::cv.glmnet(x = as.matrix(this_x_df), y = this_y, nfolds = 10,
                                            parallel = TRUE)  
      these_lasso_coefs <- coef(lassos[[i]][[j]], s = lassos[[i]][[j]]$lambda.min)
      lasso_coefs[[i]][[j]] <- as_tibble(as.matrix(these_lasso_coefs)) |> 
        mutate(Timepoint = i, Varset = j, Variable = rownames(these_lasso_coefs)) |> 
        rename(Estimate = s1)
    }
    saveRDS(glms, file = paste0(results_dir, "glms.rds"))
    saveRDS(glm_coefs, file = paste0(results_dir, "glm_coefs.rds"))
    saveRDS(lassos, file = paste0(results_dir, "lassos.rds"))
    saveRDS(lasso_coefs, file = paste0(results_dir, "lasso_coefs.rds"))
  }
}
end <- Sys.time()
cat("Elapsed time: ", format(end - start), "\n")

# get coefficients from glms, best lasso for each


cat("Data analysis (lasso, glm only) complete!\n")
