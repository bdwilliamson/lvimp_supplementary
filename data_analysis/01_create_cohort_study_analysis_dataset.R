# create analysis dataset for IMATS VIM analyses based on SRS3,
# by mimicking a prospective cohort study
library("here")
library("dplyr")
library("janitor")
library("glmnet")
library("vimp") # for stratified CV

# read in the dataset ----------------------------------------------------------
# this requires a big memory VM (*not* personal VM. Venti works, Grande ??)
data_dir <- "G:/CTRHS/IMATS/Data/SRS3 IMATS data/"
results_dir <- "G:/CTRHS/IMATS/Brian/longitudinal_vim/results/data_analysis/"
load(paste0(data_dir, "imats_srs3.Rdata"))
clean_srs3 <- imats_srs3 %>%
  clean_names()

# subset to prospective cohort study --------------------------------------
# for this approach, we'll select the closest visit to our outcome measurement time
# (90 days, 180, 270, etc.) within a 30-day window, i.e., we'll look for visits
# in the 60--120 day range, etc.
# note: would also want baseline information from everyone?
# easiest to do this by looping through and rbind'ing the dataset
srs3_18month_windows_trial <- NULL
outcome_ascertainment_times <- 1:6 * 90
outcome_window <- 30
for (i in seq_len(length(outcome_ascertainment_times))) {
  this_timepoint <- outcome_ascertainment_times[i]
  closest_date <- clean_srs3 %>% 
    mutate(days_to_time = days_since_visit1 - this_timepoint) %>% 
    filter((-1) * outcome_window <= days_to_time & days_to_time <= outcome_window) %>% 
    group_by(study_id) %>% 
    filter(abs(days_to_time) == min(abs(days_to_time))) %>% 
    mutate(timepoint = this_timepoint)
  # select those people with > 1 visit on the closest date; if none, we're done
  closest_date_with_num_obs <- closest_date %>% 
    add_count(study_id, name = "num_obs")
  multiple_visits <- closest_date_with_num_obs %>% 
    filter(num_obs > 1) %>% 
    select(-num_obs)
  single_visit <- closest_date_with_num_obs %>% 
    filter(num_obs == 1) %>% 
    select(-num_obs)
  if (nrow(multiple_visits) > 0) {
    # order to select: first, those visits marked as mental health specialty (select first)
    # second, those marked as PC1 (select first)
    # third, those marked as PC2 (select first)
    first_visit_by_type <- multiple_visits %>% 
      group_by(study_id, visit_type) %>% 
      slice(1)
    mh_visits <- first_visit_by_type %>% 
      filter(grepl("MH", visit_type)) %>% 
      ungroup()
    pc1_visits <- first_visit_by_type %>% 
      filter(grepl("PC1", visit_type), !(study_id %in% mh_visits$study_id)) %>% 
      ungroup()
    pc2_visits <- first_visit_by_type %>% 
      filter(grepl("PC2", visit_type), !(study_id %in% mh_visits$study_id), !(study_id %in% pc1_visits$study_id)) %>% 
      ungroup()
    these_data <- dplyr::bind_rows(single_visit, mh_visits, pc1_visits, pc2_visits)
  } else {
    these_data <- closest_date %>% 
      mutate(timepoint = this_timepoint)
  }
  srs3_18month_windows_trial <- dplyr::bind_rows(srs3_18month_windows_trial, these_data)
}

# save the overall dataset ------------------------------------------------
saveRDS(srs3_18month_windows_trial, file = paste0(data_dir, "imats_srs3_cohort_study_initial_dataset.rds"))

# using the cohort study dataset, subset to the variables that I need for VIM analysis ---------
srs3_18month_windows_trial <- readRDS(paste0(data_dir, "imats_srs3_cohort_study_initial_dataset.rds"))

# remove people with diagnosis of schizophrenia or bipolar disorder
srs3_trial_subset <- srs3_18month_windows_trial %>% 
  group_by(study_id) %>% 
  filter(!any(bip_visit == 1), !any(sch_visit == 1)) %>% 
  ungroup() %>% 
  select(-starts_with("bip"), -starts_with("sch"))
nrow(srs3_trial_subset)
length(unique(srs3_trial_subset$study_id))

# remove columns with variance < 0.05, and keep certain columns
vars <- srs3_trial_subset %>% 
  summarize(across(where(is.numeric), var))
var_thresh <- 0.5
sufficient_variability_vars <- names(vars)[vars > var_thresh | is.na(vars)]
srs3_trial_data_init <- srs3_trial_subset %>% 
  select(
    study_id, visit_n, visit_type, days_since_visit1, mths_since_prev, timepoint, 
    all_of(c("site", "age", "sex", "race", "hispanic")),
    charlson_score, enr_calc, census_flag, hhld_inc, coll_deg, starts_with("ins_"),
    # select diagnosis, utilization variables summarizing the total count of days with X in specific time periods
    matches("(dep|anx|bip|sch|opd|dem|add|asd|per|alc|dru|pts|eat|tbi|con|dia|ast|pai|mhi|mhe|mho|asa|lsa|osa|aip)_visit"),
    matches("(dep|anx|bip|sch|opd|dem|add|asd|per|alc|dru|pts|eat|tbi|con|dia|ast|pai|mhi|mhe|mho|asa|lsa|osa|aip)_tot_days_last[0-9]{2}m"),
    matches("(dep|anx|bip|sch|opd|dem|add|asd|per|alc|dru|pts|eat|tbi|con|dia|ast|pai|mhi|mhe|mho|asa|lsa|osa|aip)_days_per_last[0-9]{2}m"),
    # select prescription fill variables summarizing total number of months with X in specific time periods
    matches("(acv|adr|adp|ben|fga|hyp|lit|sga)_num_mths_last[0-9]{2}m"),
    matches("(acv|adr|adp|ben|fga|hyp|lit|sga)_mths_per_last[0-9]{2}m"),
    # select PHQ-8 and PHQ-9 at visit
    phq8_visit, item9_visit, 
    # select the outcome
    event90
  ) %>% 
  select(-contains("last24m"), -contains("last60m"))

srs3_trial_data_init2 <- srs3_trial_subset %>% 
  select(
    study_id, visit_n, visit_type, days_since_visit1, mths_since_prev, timepoint, 
    all_of(c("site", "age", "sex", "race", "hispanic")),
    charlson_score, enr_calc, census_flag, hhld_inc, coll_deg, starts_with("ins_"),
    all_of(sufficient_variability_vars),
    starts_with("anx"), starts_with("dep"), # anxiety, depression diagnoses
    starts_with("asa"), starts_with("lsa"), starts_with("osa"), # any, lacerative, other suicide attempt
    starts_with("aip"), event90, # accidental injury or poisoning
    # remove variables >= 2 years
    -contains("last24m", ignore.case = TRUE), -contains("last60m", ignore.case = TRUE),
    # remove delivery-related variables
    -starts_with("del_"),
    # remove 2nd way of coding information
    -ends_with("_2"),
    # remove other outcome-related information
    -starts_with("all_sui"), -starts_with("def_sui"),
    -starts_with("cen_att"), -starts_with("cen_dth"),
    -starts_with("def_death"), -event30, -event180, -contains("death"),
    # remove other PHQ information
    -(starts_with("item9") & !contains("_visit")), # remove all PHQ-9 info besides score measured at the visit
    -(starts_with("phq8") & !contains("_visit")),
    -ends_with("_2"), # remove the 2nd way of coding information
    -enr_pre_mths, -days_to_time
  )
    
ncol(srs3_trial_data_init)
names(srs3_trial_data_init)
ncol(srs3_trial_data_init2)
names(srs3_trial_data_init2)

n_events_each_timepoint <- srs3_trial_data_init %>% 
  mutate(complete = complete.cases(srs3_trial_data_init)) %>% 
  group_by(timepoint) %>% 
  summarize(n_events = sum(event90, na.rm = TRUE), 
            n_complete = sum(complete), 
            n_missing = sum(!complete),
            n_tot = n()) %>% 
  mutate(perc_missing = n_missing / n_tot * 100)
n_events_each_timepoint
readr::write_csv(n_events_each_timepoint, file = paste0(results_dir, "n_events_each_timepoint.csv"))
# timepoint n_events n_complete n_missing n_tot perc_missing
# <dbl>    <dbl>      <int>     <int> <int>        <dbl>
#   1        90      398      87527      5374 92901         5.78
# 2       180      297      64121      3977 68098         5.84
# 3       270      240      52258      3366 55624         6.05
# 4       360      200      52799      2834 55633         5.09
# 5       450      189      45958      2399 48357         4.96
# 6       540      170      41287      2203 43490         5.07
# recode character variables, add indicators for "missing" and mutate to 0 if missing, measurement otherwise
missing_cols <- is.na(srs3_trial_data_init)
which(colSums(missing_cols) > 0)
sum(is.na(srs3_trial_data_init$event90))

srs3_trial_data_init2 %>% 
  mutate(complete = complete.cases(srs3_trial_data_init2)) %>% 
  group_by(timepoint) %>% 
  summarize(n_events = sum(event90, na.rm = TRUE), 
            n_complete = sum(complete), 
            n_missing = sum(!complete),
            n_tot = n()) %>% 
  mutate(perc_missing = n_missing / n_tot * 100)
missing_cols2 <- is.na(srs3_trial_data_init2)
which(colSums(missing_cols2) > 0)
sum(is.na(srs3_trial_data_init2$event90))

srs3_trial_data <- srs3_trial_data_init %>% 
  mutate(site_hpi = as.numeric(site == "07"), # recode character vars
         visit_pc1 = as.numeric(visit_type == "PC1"),
         visit_pc2 = as.numeric(visit_type == "PC2"),
         sex_f = as.numeric(sex == "F"),
         sex_o = as.numeric(sex == "O"),
         sex_u = as.numeric(sex == "U"),
         race_as = as.numeric(race == "AS"),
         race_ba = as.numeric(race == "BA"),
         race_hp = as.numeric(race == "HP"),
         race_in = as.numeric(race == "IN"),
         race_mu = as.numeric(race == "MU"),
         race_ot = as.numeric(race == "OT"),
         race_un = as.numeric(race == "UN"), # ref: white
         hispanic = as.numeric(hispanic == "Y"),
         income_unknown = as.numeric(hhld_inc == "UNKN"),
         income_lt25 = as.numeric(hhld_inc == "LT25"),
         income_lt40 = as.numeric(hhld_inc == "LT40"), # ref: GE40
         coll_deg_unknown = as.numeric(coll_deg == "UNKN"),
         coll_deg_lt25 = as.numeric(coll_deg == "LT25"), # ref: GE25
         ins_aca_miss = as.numeric(ins_aca == "U"),
         ins_aca_num = case_when(
           ins_aca == "U" ~ 0, ins_aca == "N" ~ 1, ins_aca == "Y" ~ 2
         ),
         ins_medicaid_miss = as.numeric(ins_medicaid == "U"),
         ins_medicaid_num = case_when(
           ins_medicaid == "U" ~ 0, ins_medicaid == "N" ~ 1, ins_medicaid == "Y" ~ 2
         ),
         ins_medicare_miss = as.numeric(ins_medicare == "U"),
         ins_medicare_num = case_when(
           ins_medicare == "U" ~ 0, ins_medicare == "N" ~ 1, ins_medicare == "Y" ~ 2
         ),
         ins_privatepay_miss = as.numeric(ins_privatepay == "U"),
         ins_privatepay_num = case_when(
           ins_privatepay == "U" ~ 0, ins_privatepay == "N" ~ 1, ins_privatepay == "Y" ~ 2
         ),
         ins_statesubsidized_miss = as.numeric(ins_statesubsidized == "U"),
         ins_statesubsidized_num = case_when(
           ins_statesubsidized == "U" ~ 0, ins_statesubsidized == "N" ~ 1, ins_statesubsidized == "Y" ~ 2
         ),
         ins_selffunded_miss = as.numeric(ins_selffunded == "U"),
         ins_selffunded_num = case_when(
           ins_selffunded == "U" ~ 0, ins_selffunded == "N" ~ 1, ins_selffunded == "Y" ~ 2
         ),
         ins_highdeductible_miss = as.numeric(ins_highdeductible == "U"),
         ins_highdeductible_num = case_when(
           ins_highdeductible == "U" ~ 0, ins_highdeductible == "N" ~ 1, ins_highdeductible == "Y" ~ 2
         ),
         ins_other_miss = as.numeric(ins_other == "U"),
         ins_other_num = case_when(
           ins_other == "U" ~ 0, ins_other == "N" ~ 1, ins_other == "Y" ~ 2
         ),
         ins_commercial_miss = as.numeric(ins_commercial == "U"),
         ins_commercial_num = case_when(
           ins_commercial == "U" ~ 0, ins_commercial == "N" ~ 1, ins_commercial == "Y" ~ 2
         )
    ) %>%
  select(-site, -sex, -race, -hhld_inc, -coll_deg, -ins_aca, -ins_medicaid,
         -ins_medicare, -ins_privatepay, -ins_statesubsidized, -ins_selffunded, 
         -ins_highdeductible, -ins_commercial, -ins_other, 
         -visit_type) %>% 
  filter(!is.na(event90))
ncol(srs3_trial_data)

srs3_trial_data2 <- srs3_trial_data_init2 %>% 
  mutate(site_hpi = as.numeric(site == "07"), # recode character vars
         visit_pc1 = as.numeric(visit_type == "PC1"),
         visit_pc2 = as.numeric(visit_type == "PC2"),
         sex_f = as.numeric(sex == "F"),
         sex_o = as.numeric(sex == "O"),
         sex_u = as.numeric(sex == "U"),
         race_as = as.numeric(race == "AS"),
         race_ba = as.numeric(race == "BA"),
         race_hp = as.numeric(race == "HP"),
         race_in = as.numeric(race == "IN"),
         race_mu = as.numeric(race == "MU"),
         race_ot = as.numeric(race == "OT"),
         race_un = as.numeric(race == "UN"), # ref: white
         hispanic = as.numeric(hispanic == "Y"),
         income_unknown = as.numeric(hhld_inc == "UNKN"),
         income_lt25 = as.numeric(hhld_inc == "LT25"),
         income_lt40 = as.numeric(hhld_inc == "LT40"), # ref: GE40
         coll_deg_unknown = as.numeric(coll_deg == "UNKN"),
         coll_deg_lt25 = as.numeric(coll_deg == "LT25"), # ref: GE25
         ins_aca_miss = as.numeric(ins_aca == "U"),
         ins_aca_num = case_when(
           ins_aca == "U" ~ 0, ins_aca == "N" ~ 1, ins_aca == "Y" ~ 2
         ),
         ins_medicaid_miss = as.numeric(ins_medicaid == "U"),
         ins_medicaid_num = case_when(
           ins_medicaid == "U" ~ 0, ins_medicaid == "N" ~ 1, ins_medicaid == "Y" ~ 2
         ),
         ins_medicare_miss = as.numeric(ins_medicare == "U"),
         ins_medicare_num = case_when(
           ins_medicare == "U" ~ 0, ins_medicare == "N" ~ 1, ins_medicare == "Y" ~ 2
         ),
         ins_privatepay_miss = as.numeric(ins_privatepay == "U"),
         ins_privatepay_num = case_when(
           ins_privatepay == "U" ~ 0, ins_privatepay == "N" ~ 1, ins_privatepay == "Y" ~ 2
         ),
         ins_statesubsidized_miss = as.numeric(ins_statesubsidized == "U"),
         ins_statesubsidized_num = case_when(
           ins_statesubsidized == "U" ~ 0, ins_statesubsidized == "N" ~ 1, ins_statesubsidized == "Y" ~ 2
         ),
         ins_selffunded_miss = as.numeric(ins_selffunded == "U"),
         ins_selffunded_num = case_when(
           ins_selffunded == "U" ~ 0, ins_selffunded == "N" ~ 1, ins_selffunded == "Y" ~ 2
         ),
         ins_highdeductible_miss = as.numeric(ins_highdeductible == "U"),
         ins_highdeductible_num = case_when(
           ins_highdeductible == "U" ~ 0, ins_highdeductible == "N" ~ 1, ins_highdeductible == "Y" ~ 2
         ),
         ins_other_miss = as.numeric(ins_other == "U"),
         ins_other_num = case_when(
           ins_other == "U" ~ 0, ins_other == "N" ~ 1, ins_other == "Y" ~ 2
         ),
         ins_commercial_miss = as.numeric(ins_commercial == "U"),
         ins_commercial_num = case_when(
           ins_commercial == "U" ~ 0, ins_commercial == "N" ~ 1, ins_commercial == "Y" ~ 2
         )
  ) %>%
  select(-site, -sex, -race, -hhld_inc, -coll_deg, -ins_aca, -ins_medicaid,
         -ins_medicare, -ins_privatepay, -ins_statesubsidized, -ins_selffunded, 
         -ins_highdeductible, -ins_commercial, -ins_other, 
         -visit_type) %>% 
  filter(!is.na(event90))
ncol(srs3_trial_data2)


saveRDS(srs3_trial_data, file = paste0(data_dir, "imats_srs3_cohort_study_analysis_dataset.rds"))

# # fit lasso to everything, see what it picks
# complete_big_data <- srs3_trial_data_init %>% 
#   filter(!is.na(event90))
# big_x <- complete_big_data %>% 
#   select(-event90) %>% 
#   as.matrix()
# item9_index <- which(grepl("item9", colnames(big_x)))
# big_pen_factor <- c(rep(1, item9_index - 1), 0, rep(1, ncol(big_x) - item9_index))
# big_lasso <- glmnet::cv.glmnet(x = big_x, y = complete_big_data$event90, 
#                                family = "binomial",
#                                penalty.factor = big_pen_factor, nfolds = 10)

# check manually-selected vars using glm
timepoints <- sort(unique(srs3_trial_data$timepoint))
K <- 10
aucs_glm <- matrix(nrow = length(timepoints), ncol = K)
aucs2 <- aucs_glm2 <- matrix(nrow = length(timepoints), ncol = K)
set.seed(20230403)
seeds <- round(runif(length(timepoints), 1e4, 1e5))
for (i in 1:length(timepoints)) {
  these_data <- srs3_trial_data %>% 
    filter(timepoint == timepoints[i], !is.na(event90)) %>% 
    select(-timepoint, -visit_n, -study_id) %>% 
    na.omit()
  set.seed(seeds[i])
  # these_folds <- sample(seq_len(K), nrow(these_data), replace = TRUE)
  these_folds <- vimp::make_folds(y = these_data$event90, V = K, stratified = TRUE)
  for (k in seq_len(K)) {
    this_train <- these_data[these_folds != k, ]
    this_test <- these_data[these_folds == k, ]
    fit <- glm("event90 ~ .", data = this_train, family = "binomial")
    preds <- predict(fit, newdata = this_test, type = "response")
    aucs_glm[i, k] <- cvAUC::cvAUC(predictions = preds, labels = this_test$event90)$fold.AUC
  }
}
rowMeans(aucs_glm)

# check manually-selected vars using lasso
timepoints <- sort(unique(srs3_trial_data$timepoint))
K <- 10
aucs_lasso <- matrix(nrow = length(timepoints), ncol = K)
lasso_k <- 5
set.seed(20230403)
seeds <- round(runif(length(timepoints), 1e4, 1e5))
for (i in 1:length(timepoints)) {
  these_data <- srs3_trial_data %>% 
    filter(timepoint == timepoints[i], !is.na(event90)) %>% 
    select(-timepoint, -visit_n, -study_id) %>% 
    na.omit()
  set.seed(seeds[i])
  # these_folds <- sample(seq_len(K), nrow(these_data), replace = TRUE)
  these_folds <- vimp::make_folds(y = these_data$event90, V = K, stratified = TRUE)
  for (k in seq_len(K)) {
    this_train <- these_data[these_folds != k, ]
    this_test <- these_data[these_folds == k, ]
    this_x <- as.matrix(this_train[, -which(names(this_train) == "event90")])
    item9_index <- which(grepl("item9", colnames(this_x)))
    penalty_factor <- c(rep(1, item9_index - 1), 0, rep(1, ncol(this_x) - item9_index))
    this_y <- this_train$event90
    this_glmnet_folds <- vimp::make_folds(y = this_y, V = lasso_k, stratified = TRUE)
    fit <- glmnet::cv.glmnet(x = this_x, y = this_y, family = "binomial",
                             penalty.factor = penalty_factor, nfolds = lasso_k, foldid = this_glmnet_folds)
    test_x <- this_test %>%
      select(-event90) %>%
      as.matrix()
    preds <- predict(fit, newx = test_x, s = fit$lambda.min, type = "response")
    aucs_lasso[i, k] <- cvAUC::cvAUC(predictions = preds, labels = this_test$event90)$fold.AUC
  }
}
rowMeans(aucs_lasso)

# fit a lasso at all timepoints; pick only those with nonzero coefs
set.seed(20230512)
K <- 10
big_x <- srs3_trial_data2 %>% 
  select(-timepoint, -visit_n, -study_id, -event90) %>% 
  as.matrix()
item9_index <- which(grepl("item9", colnames(big_x)))
big_pen_factor <- c(rep(1, item9_index - 1), 0, rep(1, ncol(big_x) - item9_index))
big_folds <- vimp::make_folds(y = srs3_trial_data2$event90, V = K, stratified = TRUE)
big_lasso <- glmnet::cv.glmnet(x = big_x, y = srs3_trial_data2$event90,
                               family = "binomial",
                               penalty.factor = big_pen_factor, nfolds = 10, foldid = big_folds)
lasso_vars <- as.logical(abs(coef(big_lasso$glmnet.fit, s = big_lasso$lambda.min)) > 0)[-1]
lasso_var_names <- colnames(big_x)[lasso_vars]
# save off a dataset with variables selected by lasso
srs3_trial_data_lasso <- srs3_trial_data2 %>% 
  select(timepoint, visit_n, study_id, event90, all_of(lasso_var_names))
saveRDS(srs3_trial_data_lasso, file = paste0(data_dir, "imats_srs3_cohort_study_analysis_dataset_lasso_selected.rds"))

# check to see which variables are in one but not the other
nms_preselected <- names(srs3_trial_data)
nms_lasso <- names(srs3_trial_data_lasso)
setdiff(nms_preselected, nms_lasso)
setdiff(nms_lasso, nms_preselected)
intersect(nms_lasso, nms_preselected)

K <- 10
# check lasso-selected vars using glm
for (i in 1:length(timepoints)) {
  these_data <- srs3_trial_data2 %>% 
    filter(timepoint == timepoints[i], !is.na(event90)) %>% 
    select(all_of(lasso_var_names), event90)
  set.seed(seeds[i])
  these_folds <- vimp::make_folds(y = these_data$event90, V = K, stratified = TRUE)
  for (k in seq_len(K)) {
    this_train <- these_data[these_folds != k, ]
    this_test <- these_data[these_folds == k, ]
    fit <- glm("event90 ~ .", data = this_train, family = "binomial")
    preds <- predict(fit, newdata = this_test, type = "response")
    aucs_glm2[i, k] <- cvAUC::cvAUC(predictions = preds, labels = this_test$event90)$fold.AUC
  }
}
rowMeans(aucs_glm2)

# look at the 500-variable dataset
timepoints <- sort(unique(srs3_trial_data$timepoint))
set.seed(20230403)
seeds <- round(runif(length(timepoints), 1e4, 1e5))
K <- 5
lasso_k <- 5
aucs <- matrix(nrow = length(timepoints), ncol = K)
for (i in 1:length(timepoints)) {
  these_data <- srs3_trial_data2 %>% 
    filter(timepoint == timepoints[i], !is.na(event90)) %>% 
    select(-timepoint, -visit_n, -study_id) %>% 
    na.omit()
  set.seed(seeds[i])
  these_folds <- vimp::make_folds(y = these_data$event90, V = K, stratified = TRUE)
  for (k in seq_len(K)) {
    this_train <- these_data[these_folds != k, ]
    this_test <- these_data[these_folds == k, ]
    this_x <- as.matrix(this_train[, -which(names(this_train) == "event90")])
    item9_index <- which(grepl("item9", colnames(this_x)))
    penalty_factor <- c(rep(1, item9_index - 1), 0, rep(1, ncol(this_x) - item9_index))
    this_y <- this_train$event90
    this_glmnet_folds <- vimp::make_folds(y = this_y, V = lasso_k, stratified = TRUE)
    fit <- glmnet::cv.glmnet(x = this_x, y = this_y, family = "binomial",
                             penalty.factor = penalty_factor, nfolds = lasso_k, foldid = this_glmnet_folds)
    test_x <- this_test %>%
      select(-event90) %>%
      as.matrix()
    preds <- predict(fit, newx = test_x, s = fit$lambda.min, type = "response")
    aucs[i, k] <- cvAUC::cvAUC(predictions = preds, labels = this_test$event90)$fold.AUC
  }
}
rowMeans(aucs)
coef(fit$glmnet.fit, s = fit$lambda.min)
