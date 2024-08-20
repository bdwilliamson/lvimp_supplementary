# plots and tables of results

# load required packages and functions -----------------------------------------
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("dplyr")
library("tidyr")
library("here")
library("rprojroot")
this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(proj_root, "/code/sims/utils.R"))
source(paste0(proj_root, "/code/data_analysis/00_utils.R"))

# read in the results ----------------------------------------------------------
data_dir <- "G:/CTRHS/IMATS/Data/SRS3 IMATS data/"
results_dir <- "G:/CTRHS/IMATS/Brian/longitudinal_vim/results/data_analysis/"

output <- readRDS(file = paste0(results_dir, "lvim_output.rds")) 
glm_coefs <- readRDS(file = paste0(results_dir, "glm_coefs.rds"))
lasso_coefs <- readRDS(file = paste0(results_dir, "lasso_coefs.rds"))

# set up the variable sets
analysis_dataset <- readRDS(paste0(data_dir, "imats_srs3_cohort_study_analysis_dataset.rds"))
y <- analysis_dataset$event90
X <- analysis_dataset %>% 
  select(-study_id, -visit_n, -timepoint, -event90)
timepoints <- unique(analysis_dataset$timepoint)
num_timepoints <- length(timepoints)
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
  list(c(varset_demographics, varset_prior_self_harm, varset_phq9)), # (6) add PHQi9 to baseline model
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson)), # (7) add charlson to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_phq9)), # (8) add PHQi9 to model (4)
  list(c(varset_demographics, varset_prior_self_harm, varset_charlson, varset_phq9)), # (9) add PHQi9 to model (5)
  list(c(varset_demographics, varset_prior_self_harm, varset_diagnoses, varset_charlson, varset_phq9)), # (10) all variables
  list(varset_phq9), # (11) PHQi9 compared to nothing, corresponds to variable set 2 in the paper
  list(varset_age), # (12) age compared to nothing
  list(varset_age_sex), # (13) age and sex compared to nothing, corresponds to variable set 3 in the paper
  list(varset_age_sex_race_ethnicity), # (14) age, sex, race and ethnicity compared to nothing
  list(c(varset_age, varset_phq9)), # (15) PHQi9 compared to age alone (12)
  list(c(varset_age_sex, varset_phq9)), # (16) PHQi9 compared to age and sex (13), corresponds to variable set 4 in the paper
  list(c(varset_age_sex_race_ethnicity, varset_phq9)), # (17) PHQi9 compared to age, sex, race and ethnicity (14)
  list(c(varset_age_sex, varset_prior_self_harm)), # (18) corresponds to variable set 5 in the paper
  list(c(varset_age_sex, varset_prior_self_harm, varset_phq9)) # (19) corresponds to variable set 6 in the paper
)
num_varsets <- length(varsets)
easy_varsets <- 1:num_varsets
paper_varsets <- c(1, 7:15, 2, 16, 3, 17, 18, 4, 19, 5, 6)
paper_varset_df <- tibble::tibble(varset_num = easy_varsets, paper_varset_num = paper_varsets)
nice_varsets <- c("No variables", 
                  "Demographic variables", 
                  "Demographic and prior self-harm variables",
                  "Demographic, prior self-harm, and diagnosis & utilization variables",
                  "Demographic and prior self-harm variables and Charlson score",
                  "Demographic and prior self-harm variables and PHQi9",
                  "Demographic, prior self-harm, and diagnosis & utilization variables, and Charlson score",
                  "Demographic, prior self-harm, and diagnosis & utilization variables, and PHQi9",
                  "Demographic and prior self-harm variables, Charlson score, and PHQi9",
                  "All variables",
                  "PHQi9",
                  "Age",
                  "Age and sex",
                  "Age, sex, and race or ethnicity",
                  "Age and PHQi9",
                  "Age, sex, and PHQi9",
                  "Age, sex, race or ethnicity, and PHQi9",
                  "Age, sex, and prior self-harm variables",
                  "Age, sex, prior self-harm variables, and PHQi9")
nice_output_init <- output %>% 
  mutate(varset_num = get_varset_num(varset)) %>% 
  mutate(nice_varset = get_nice_varset(varset_num, nice_varsets),
         algo_fct = factor(algo, labels = c("MEAN", "GLM", "LASSO", "RF", "Discrete SL", "SL", "XGB"),
                           levels = c("mean", "glm", "glmnet", "rf", "Discrete SL", "SL", "xgb"),
                           ordered = TRUE))
nice_output <- nice_output_init |> 
  left_join(paper_varset_df, by = "varset_num") |> 
  mutate(varset_fct = factor(paper_varset_num, levels = paper_varsets, labels = nice_varsets, ordered = TRUE)) %>% 
  select(-varset) %>% 
  filter(!grepl("mean", algo))

# create plots and tables ------------------------------------------------------
dodge_width <- 0.875
legend_text_size <- 10
axis_text_size <- 10
title_text_size <- 12
fig_width <- 9
fig_height <- 6
point_size <- 1.5
# plot predictiveness of variable sets 1, 2, 3 (original numbering)
# or variable sets 1, 7, 8 (paper numbering)
varset_178_pred_plot <- plot_trajectory(output = nice_output, varset = c(1, 7, 8),
                                         est_type = "predictiveness")
# plot predictiveness of variable set 4
# variable set 9 in paper ordering
varset_9_pred_plot <- plot_trajectory(output = nice_output, varset = 9,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 5
# variable set 10 in paper ordering
varset_10_pred_plot <- plot_trajectory(output = nice_output, varset = 10,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 6
# variable set 11 in paper ordering
varset_11_pred_plot <- plot_trajectory(output = nice_output, varset = 11,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 7
# variable set 12 in paper ordering
varset_12_pred_plot <- plot_trajectory(output = nice_output, varset = 12,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 8
# variable set 13 in paper ordering
varset_13_pred_plot <- plot_trajectory(output = nice_output, varset = 13,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 9
# variable set 14 in paper ordering
varset_14_pred_plot <- plot_trajectory(output = nice_output, varset = 14,
                                      est_type = "predictiveness")
# plot predictiveness of variable set 10
# variable set 15 in paper ordering
varset_15_pred_plot <- plot_trajectory(output = nice_output, varset = 15,
                                      est_type = "predictiveness")

# plot predictiveness of variable sets 11, 12, 13, 14 (each is compared to null set)
# variable sets 2, 16, 3, 17 in paper ordering
varset_2_3_16_17_pred_plot <- plot_trajectory(output = nice_output, varset = c(2, 3, 16, 17),
                                             est_type = "predictiveness")
# variable sets 18, 4, 19 in paper ordering
varset_4_18_19_pred_plot <- plot_trajectory(output = nice_output, varset = c(4, 18, 19),
                                             est_type = "predictiveness")

# plot VIM of variable sets 
varset_2_vim_plot <- plot_trajectory(output = nice_output, varset = 2,
                                      est_type = "vim")
varset_3_vim_plot <- plot_trajectory(output = nice_output, varset = 3,
                                       est_type = "vim")
varset_4_vim_plot <- plot_trajectory(output = nice_output, varset = 4,
                                     est_type = "vim")
varset_5_vim_plot <- plot_trajectory(output = nice_output, varset = 5,
                                     est_type = "vim")
varset_6_vim_plot <- plot_trajectory(output = nice_output, varset = 6,
                                     est_type = "vim")
varset_7_vim_plot <- plot_trajectory(output = nice_output, varset = 7,
                                     est_type = "vim")
varset_8_vim_plot <- plot_trajectory(output = nice_output, varset = 8,
                                     est_type = "vim")
varset_9_vim_plot <- plot_trajectory(output = nice_output, varset = 9,
                                     est_type = "vim")
varset_10_vim_plot <- plot_trajectory(output = nice_output, varset = 10,
                                     est_type = "vim")
varset_11_vim_plot <- plot_trajectory(output = nice_output, varset = 11,
                                      est_type = "vim")
varset_12_vim_plot <- plot_trajectory(output = nice_output, varset = 12,
                                      est_type = "vim")
varset_13_vim_plot <- plot_trajectory(output = nice_output, varset = 13,
                                      est_type = "vim")
varset_14_vim_plot <- plot_trajectory(output = nice_output, varset = 14,
                                      est_type = "vim")
varset_15_vim_plot <- plot_trajectory(output = nice_output, varset = 15,
                                      est_type = "vim")
varset_16_vim_plot <- plot_trajectory(output = nice_output, varset = 16,
                                      est_type = "vim")
varset_17_vim_plot <- plot_trajectory(output = nice_output, varset = 17,
                                      est_type = "vim")
varset_18_vim_plot <- plot_trajectory(output = nice_output, varset = 18,
                                      est_type = "vim")
varset_19_vim_plot <- plot_trajectory(output = nice_output, varset = 19,
                                      est_type = "vim")

# tables with predictiveness or VIM for each timepoint
all_pred_table <- create_overall_table(output = nice_output, varset = 1:19,
                                       est_type = "predictiveness") 
all_vim_table <- create_overall_table(output = nice_output, varset = 1:19,
                                      est_type = "vim")
vim_lag1_table <- all_vim_table %>% 
  mutate(var = SE ^ 2) %>% 
  select(-`95% CI`, -SE, -`p-value`) %>% 
  pivot_wider(names_from = Timepoint, values_from = c(`Point Estimate`, var)) %>% 
  mutate(vim_est_lag_2 = `Point Estimate_2` - `Point Estimate_1`,
         vim_var_lag_2 = var_1 + var_2,
         vim_est_lag_3 = `Point Estimate_3` - `Point Estimate_2`,
         vim_var_lag_3 = var_2 + var_3,
         vim_est_lag_4 = `Point Estimate_4` - `Point Estimate_3`,
         vim_var_lag_4 = var_3 + var_4,
         vim_est_lag_5 = `Point Estimate_5` - `Point Estimate_4`,
         vim_var_lag_5 = var_4 + var_5,
         vim_est_lag_6 = `Point Estimate_6` - `Point Estimate_5`,
         vim_var_lag_6 = var_5 + var_6) %>% 
  select(-starts_with("Point Estimate"), -starts_with("var_")) %>% 
  pivot_longer(cols = starts_with("vim"), 
               names_to = c("type", "Timepoint"),
               names_pattern = "(.*)_lag_(.*)") %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(se = sqrt(vim_var), cil = vim_est - qnorm(0.975) * se,
         ciu = vim_est + qnorm(0.975) * se,
         `95% CI` = paste0("[", round(cil, 3), ", ", round(ciu, 3), "]"),
         `p-value` = 1 - pnorm(abs(vim_est / se))) %>% 
  select(-cil, -ciu) %>% 
  rename(`Point Estimate` = vim_est, SE = se) %>% 
  mutate(Timepoint = paste0(Timepoint, " - ", as.numeric(Timepoint) - 1))
  

# table with average, slope (and intercept) of line
all_pred_summary_table <- create_overall_table(output = nice_output, varset = 1:19,
                                               est_type = "predictiveness",
                                               trajectory = FALSE) |> 
  filter(!grepl("AUTC", `Summary measure`))
all_vim_summary_table <- create_overall_table(output = nice_output, varset = 1:19,
                                              est_type = "vim",
                                              trajectory = FALSE) |> 
  filter(!grepl("AUTC", `Summary measure`))

# write out plots
pred_plot_nms <- c("178", as.character(9:15), "2_3_16_17", "4_18_19")
for (nm in pred_plot_nms) {
  eval(parse(text = paste0("ggsave(filename = paste0(results_dir, 'pred_trajectory_", nm, ".png'), ", 
                           "plot = varset_", nm, "_pred_plot, width = 8, height = 11, units = 'in', dpi = 300)")))
}
vim_plot_nms <- as.character(2:num_varsets)
for (nm in vim_plot_nms) {
  eval(parse(text = paste0("ggsave(filename = paste0(results_dir, 'vim_trajectory_", nm, ".png'), ", 
                           "plot = varset_", nm, "_vim_plot, width = 8, height = 11, units = 'in', dpi = 300)")))
}

# write out tables
readr::write_csv(all_pred_table, file = paste0(results_dir, "pred_trajectory.csv"))
readr::write_csv(all_vim_table, file = paste0(results_dir, "vim_trajectory.csv"))
readr::write_csv(vim_lag1_table, file = paste0(results_dir, "vim_lag1.csv"))

readr::write_csv(all_pred_summary_table, file = paste0(results_dir, "pred_summary.csv"))
readr::write_csv(all_vim_summary_table, file = paste0(results_dir, "vim_summary.csv"))

# main figure: showing VIM over time for several comparisons
varsets_of_interest <- c(2:6)
# varsets_of_interest <- c(11, 15, 16, 17, 6, 10)
# other_varsets_of_interest <- c(2, 18, 4, 19, 11, 15)
plot_tib <- nice_output |> 
  filter(paper_varset_num %in% varsets_of_interest, 
         designation %in% paste0("vim-", 1:num_timepoints)) |>
  mutate(timepoint = gsub("vim-", "", designation),
         varset_num_fct = factor(paper_varset_num, levels = varsets_of_interest,
                                 labels = c("Add-in: variable set 2 vs 1", 
                                            "Add-in: variable set 3 vs 1",
                                            "Add-in: variable set 4 vs 3",
                                            "Add-in: variable set 5 vs 3",
                                            "Add-in: variable set 6 vs 5"), ordered = TRUE),
         truncated_est = pmax(est, 0))
  
main_fig_all_estimators <- plot_tib |> 
  ggplot(aes(x = timepoint, y = truncated_est)) +
  geom_point(position = position_dodge(dodge_width)) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(dodge_width)) +
  scale_color_viridis_d(begin = 0, end = 0.75) +
  labs(x = "Timepoint", y = "Estimated VIM",
       color = "Estimator", shape = "Estimator") +
  ylim(c(-0.01, 0.21)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(cols = vars(varset_num_fct), rows = vars(algo_fct)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        strip.text = element_text(size = 8))
ggsave(filename = paste0(results_dir, "vims_of_interest_over_time_all_estimators.png"),
       main_fig_all_estimators, width = 8.5, height = 6, units = "in",
       dpi = 300)

main_fig <- plot_tib |> 
  filter(algo %in% c("glm", "SL")) |> 
  ggplot(aes(x = timepoint, y = truncated_est)) +
  geom_point(position = position_dodge(dodge_width)) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(dodge_width)) +
  scale_color_viridis_d(begin = 0, end = 0.75) +
  labs(x = "Timepoint", y = "Estimated VIM",
       color = "Estimator", shape = "Estimator") +
  ylim(c(-0.01, 0.21)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(cols = vars(varset_num_fct), rows = vars(algo_fct)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        strip.text = element_text(size = 8))
ggsave(filename = paste0(results_dir, "vims_of_interest_over_time.png"),
       main_fig, width = 8.5, height = 3, units = "in",
       dpi = 300)

plot_tib_presentation <- plot_tib |> 
  filter(algo %in% c("SL")) |> 
  mutate(nicer_varset = factor(
    case_when(
      varset_fct == "PHQi9" ~ "PHQi9 vs no variables",
      varset_fct == "Age and sex" ~ "Age and sex vs no variables",
      varset_fct == "Age, sex, and PHQi9" ~ "PHQi9 vs age and sex",
      varset_fct == "Age, sex, and prior self-harm variables" ~ "Prior self-harm vs age and sex",
      varset_fct == "Age, sex, prior self-harm variables, and PHQi9" ~ "PHQi9 vs age, sex, and prior self-harm"
    ), labels = c("Age and sex vs no variables", "PHQi9 vs no variables", "Prior self-harm vs age and sex",
                  "PHQi9 vs age and sex", "PHQi9 vs age, sex, and prior self-harm"), 
    levels = c("Age and sex vs no variables", "PHQi9 vs no variables", "Prior self-harm vs age and sex",
               "PHQi9 vs age and sex", "PHQi9 vs age, sex, and prior self-harm"), ordered = TRUE
  ))
main_presentation_fig <- plot_tib_presentation |> 
  ggplot(aes(x = timepoint, y = truncated_est)) +
  geom_point(position = position_dodge(dodge_width)) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(dodge_width)) +
  scale_color_viridis_d(begin = 0, end = 0.75) +
  labs(x = "Timepoint", y = "Estimated VIM",
       color = "Estimator", shape = "Estimator") +
  ylim(c(-0.01, 0.21)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(vars(nicer_varset), nrow = 2, ncol = 3) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        strip.text = element_text(size = 8))
ggsave(filename = paste0(results_dir, "vims_of_interest_over_time_presentation.png"),
       main_presentation_fig, width = 8.5, height = 4.5, units = "in",
       dpi = 300)

main_tab <- nice_output |> 
  filter(paper_varset_num %in% varsets_of_interest, designation %in% c("vim-average", "vim-trend-slope")) |> 
  mutate(designation = gsub("vim-", "", designation),
         varset_num_fct = factor(paper_varset_num, levels = varsets_of_interest,
                                 labels = c("PHQi9 versus no variables", 
                                            "Age and sex versus no variables",
                                            "PHQi9 versus age and sex",
                                            "Prior self-harm variables versus age and sex",
                                            "PHQi9 versus age, sex, and prior self-harm variables"), ordered = TRUE)) |> 
  arrange(algo_fct, designation, varset_num_fct) |> 
  select(designation, varset_num_fct, algo_fct, est, se, cil, ciu, p_value) |> 
  mutate(est = sprintf("%.3f", ifelse(abs(est) < 0.001, 0, est)), se = sprintf("%.3f", se), 
         ci = paste0("[", sprintf("%.3f", cil), ", ", sprintf("%.3f", ciu), "]"),
         p_value = ifelse(p_value < 1e-3, "< 0.001", sprintf("%.3f", p_value))) |> 
  rename(Estimator = algo_fct, `VIM type: comparison` = varset_num_fct,
         `Summary` = designation, 
         `Estimate` = est, 
         `SE` = se,
         `95\\% CI` = ci, `p-value` = p_value) |> 
  select(Summary, `VIM type: comparison`, `Estimator`, `Estimate`, `SE`, `95\\% CI`, `p-value`)

main_tab_presentation <- main_tab |> 
  filter(Estimator == "SL") |> 
  select(-Estimator) |> 
  rename(`95% CI` = `95\\% CI`) |> 
  mutate(Summary = case_when(
    Summary == "average" ~ "Mean",
    Summary == "trend-slope" ~ "Trend: slope"
  ), Comparison = factor(`VIM type: comparison`, 
                         labels = c("Age and sex versus no variables", 
                                    "PHQi9 versus no variables",
                                    "Prior self-harm variables versus age and sex",
                                    "PHQi9 versus age and sex",
                                    "PHQi9 versus age, sex, and prior self-harm variables"), 
                         levels = c("Age and sex versus no variables", 
                                    "PHQi9 versus no variables",
                                    "Prior self-harm variables versus age and sex",
                                    "PHQi9 versus age and sex",
                                    "PHQi9 versus age, sex, and prior self-harm variables"), ordered = TRUE)) |> 
  arrange(Summary, Comparison) |> 
  select(-`VIM type: comparison`) |> 
  select(Summary, Comparison, everything())
saveRDS(main_tab_presentation, file = paste0(results_dir, "vim_summaries_presentation.rds"))
# other variable sets
# 5 -> 10, 7 -> 12, 8 -> 13, 9 -> 14
# c(2, 18, 4, 19, 11, 15)
other_varsets <- c(10:15, 18:19)
supp_tab <- nice_output |> 
  filter(paper_varset_num %in% other_varsets, designation %in% c("vim-average", "vim-trend-slope")) |> 
  mutate(designation = gsub("vim-", "", designation),
         varset_num_fct = factor(paper_varset_num, levels = other_varsets,
                                 labels = c("Charlson vs demographic, prior self-harm variables", 
                                            "PHQi9 vs demographic, prior self-harm variables",
                                            "Charlson vs demographic, prior self-harm, diagnosis \\& utilization variables",
                                            "PHQi9 vs demographic, prior self-harm, diagnosis \\& utilization variables",
                                            "PHQi9 vs demographic, prior self-harm variables and Charlson score",
                                            "PHQi9 vs all other variables",
                                            "PHQi9 vs age",
                                            "PHQi9 vs age, sex, and race or ethnicity"), ordered = TRUE)) |> 
  arrange(algo_fct, designation, varset_num_fct) |> 
  select(designation, varset_num_fct, algo_fct, est, se, cil, ciu, p_value) |> 
  mutate(est = sprintf("%.3f", ifelse(abs(est) < 0.001, 0, est)), se = sprintf("%.3f", se), 
         ci = paste0("[", sprintf("%.3f", cil), ", ", sprintf("%.3f", ciu), "]"),
         p_value = ifelse(p_value < 1e-3, "< 0.001", sprintf("%.3f", p_value))) |> 
  rename(Estimator = algo_fct, `VIM type: comparison` = varset_num_fct,
         `Summary` = designation, 
         `Estimate` = est, 
         `SE` = se,
         `95\\% CI` = ci, `p-value` = p_value) |> 
  select(Summary, `VIM type: comparison`, `Estimator`, `Estimate`, `SE`, `95\\% CI`, `p-value`)


font_size <- 9
main_tab |> 
  filter(Estimator == "GLM" | Estimator == "SL") |> 
  select(-Summary, -Estimator) |> 
  knitr::kable(digits = 3, format = "latex", booktabs = TRUE, escape = FALSE,
               caption = paste0("Estimates of the average VIM and slope of the linear trend in VIM over the time series, ",
                                " considering the importance of PHQi9 when compared to other variables.",
                                " Comparisons are: PHQi9 vs no variables (comparing variable set 2 to 1),",
                                " age and sex vs no variables (comparing variable set 3 to 1),",
                                " PHQi9 vs age and sex (comparing variable set 4 to 3),",
                                " prior self-harm variables vs age and sex (comparing variable set 5 to 3),",
                                " PHQi9 vs age, sex, and prior self-harm variables (comparing variable set 6 to 5).",
                                " Estimates are shown for logistic regression (GLM) and Super Learner (SL).",
                                " \\label{tab:data_analysis_vim_summaries}"),
               linesep = "") |> 
  kableExtra::kable_styling(font_size = font_size, 
                            latex_options = c("scale_down")) |> 
  kableExtra::pack_rows("Average (GLM)", 1, 5) |> 
  kableExtra::pack_rows("Trend - slope (GLM)", 6, 10) |> 
  kableExtra::pack_rows("Average (SL)", 11, 15) |> 
  kableExtra::pack_rows("Trend - slope (SL)", 16, 20) |> 
  kableExtra::save_kable(file = paste0(results_dir, "data_analysis_vim_summaries.tex"))

main_tab |> 
  select(-Summary, -Estimator) |> 
  knitr::kable(digits = 3, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
               caption = paste0("Estimates of the average VIM and slope of the linear trend in VIM over the time series, ",
                                " considering the importance of PHQi9 when compared to other variables.",
                                " Comparisons are: PHQi9 vs no variables (comparing variable set 2 to 1),",
                                " age and sex vs no variables (comparing variable set 3 to 1),",
                                " PHQi9 vs age and sex (comparing variable set 4 to 3),",
                                " prior self-harm variables vs age and sex (comparing variable set 5 to 3),",
                                " PHQi9 vs age, sex, and prior self-harm variables (comparing variable set 6 to 5).",
                                " Estimates are shown for logistic regression (GLM), lasso (LASSO),",
                                " random forests (RF), the discrete Super Learner (Discrete SL), Super Learner (SL),",
                                " and boosted trees (XGB).",
                                " \\label{tab:data_analysis_vim_summaries_all}"),
               linesep = "") |> 
  kableExtra::kable_styling(font_size = font_size, 
                            latex_options = c("repeat_header"),
                            repeat_header_text = "Estimates of the average VIM and slope of the linear trend in VIM over the time series, considering the importance of PHQi9 when compared to other variables. \\textit{(continued)}",
                            repeat_header_method = "replace") |> 
  kableExtra::pack_rows("Average (GLM)", 1, 5) |> 
  kableExtra::pack_rows("Trend - slope (GLM)", 6, 10) |> 
  kableExtra::pack_rows("Average (LASSO)", 11, 15) |> 
  kableExtra::pack_rows("Trend - slope (LASSO)", 16, 20) |> 
  kableExtra::pack_rows("Average (RF)", 21, 25) |> 
  kableExtra::pack_rows("Trend - slope (RF)", 26, 30) |> 
  kableExtra::pack_rows("Average (Discrete SL)", 31, 35) |> 
  kableExtra::pack_rows("Trend - slope (Discrete SL)", 36, 40) |> 
  kableExtra::pack_rows("Average (SL)", 41, 45) |> 
  kableExtra::pack_rows("Trend - slope (SL)", 46, 50) |> 
  kableExtra::pack_rows("Average (XGB)", 51, 55) |> 
  kableExtra::pack_rows("Trend - slope (XGB)", 56, 60) |> 
  kableExtra::save_kable(file = paste0(results_dir, "data_analysis_vim_summaries_all.tex"))

supp_tab |> 
  select(-Summary, -Estimator) |> 
  knitr::kable(digits = 3, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
               caption = paste0("Estimates of the average VIM and slope of the linear trend in VIM over the time series, ",
                                " considering the importance of PHQi9 or Charlson score when compared to other variables.",
                                " Comparisons are:", 
                                " Charlson vs demographic and prior self-harm variables (comparing variable set 10 to 8),",
                                " PHQi9 vs demographic and prior self-harm variables (comparing variable set 11 to 8),",
                                " Charlson vs demographic, prior self-harm, and diagnosis \\& utilization variables (comparing variable set 12 to 9),",
                                " PHQi9 vs demographic, prior self-harm, and diagnosis \\& utilization variables (comparing variable set 13 to 9),",
                                " PHQi9 vs demographic, prior self-harm variables and Charlson score (comparing variable set 14 to 10),",
                                " PHQi9 vs all other variables (comparing variable set 15 to 12),",
                                " PHQi9 vs age (comparing variable set 18 to 16),",
                                " PHQi9 vs age, sex, and race or ethnicity (comparing variable set 19 to 17).",
                                " Estimates are shown for logistic regression (GLM), lasso (LASSO),",
                                " random forests (RF), the discrete Super Learner (Discrete SL), Super Learner (SL),",
                                " and boosted trees (XGB).",
                                " \\label{tab:data_analysis_vim_summaries_all_supp}"),
               linesep = "") |> 
  kableExtra::kable_styling(font_size = font_size, 
                            latex_options = c("repeat_header"),
                            repeat_header_text = paste0("Estimates of the average VIM and slope of the linear trend in VIM over the time series, ",
                                                        " considering the importance of PHQi9 or Charlson score when compared to other variables. \\textit{(continued)}"),
                            repeat_header_method = "replace") |> 
  kableExtra::pack_rows("Average (GLM)", 1, 8) |> 
  kableExtra::pack_rows("Trend - slope (GLM)", 9, 16) |> 
  kableExtra::pack_rows("Average (LASSO)", 17, 24) |> 
  kableExtra::pack_rows("Trend - slope (LASSO)", 25, 32) |> 
  kableExtra::pack_rows("Average (RF)", 33, 40) |> 
  kableExtra::pack_rows("Trend - slope (RF)", 41, 48) |> 
  kableExtra::pack_rows("Average (Discrete SL)", 49, 56) |> 
  kableExtra::pack_rows("Trend - slope (Discrete SL)", 57, 64) |> 
  kableExtra::pack_rows("Average (SL)", 65, 72) |> 
  kableExtra::pack_rows("Trend - slope (SL)", 73, 80) |> 
  kableExtra::pack_rows("Average (XGB)", 81, 88) |> 
  kableExtra::pack_rows("Trend - slope (XGB)", 89, 96) |> 
  kableExtra::save_kable(file = paste0(results_dir, "data_analysis_vim_summaries_all_supp.tex"))
  

# investigate glm, lasso for variable set with PHQ and baseline variables -----------------
glm_coef_df <- do.call(rbind, do.call(Map, c(f = rbind, glm_coefs)))
nice_glm_coefs <- glm_coef_df |> 
  mutate(nice_varset = get_nice_varset(Varset, nice_varsets),
         varset_fct = factor(Varset, levels = 1:num_varsets, labels = nice_varsets, ordered = TRUE))
# drop the intercept; look at coefficients over time for the following models:
# (11) PHQi9 alone
# (15) PHQi9 and age
# (16) PHQi9, age, sex
# (17) PHQi9, age, sex, race or ethnicity
# (6)  PHQi9 + baseline variables
nice_glm_coefs |> 
  filter(Varset %in% c(6, 11, 15, 16, 17)) |> 
  group_by(Varset, Variable) |> 
  summarize(across(Estimate, .fns = list(mn = ~ mean(.x), sd = ~ sd(.x), 
                                      min = ~ min(.x), max = ~ max(.x),
                                      pct25 = ~ quantile(.x, .25),
                                      pct50 = ~ median(.x),
                                      pct75 = ~ quantile(.x, .75)),
                .names = "{.fn}_{.col}")) |> 
  readr::write_csv(file = paste0(results_dir, "summ_glm_coefs.csv"))

phq9_glm_plot <- nice_glm_coefs |> 
  filter(Varset %in% c(6, 11, 15, 16, 17), grepl("item9", Variable)) |> 
  ggplot(aes(x = Timepoint, y = Estimate, color = nice_varset, shape = nice_varset)) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(color = "Variable group", shape = "Variable group") 
ggsave(filename = paste0(results_dir, "plot_glm_phq9.png"), width = 10, height = 5,
       plot = phq9_glm_plot)

lasso_coef_df <- data.table::rbindlist(lapply(lasso_coefs, data.table::rbindlist, fill = TRUE))
nice_lasso_coefs <- lasso_coef_df |> 
  mutate(nice_varset = get_nice_varset(Varset, nice_varsets),
         varset_fct = factor(Varset, levels = 1:num_varsets, labels = nice_varsets, ordered = TRUE))
# drop the intercept; look at coefficients over time for the following models:
# (11) PHQi9 alone
# (15) PHQi9 and age
# (16) PHQi9, age, sex
# (17) PHQi9, age, sex, race or ethnicity
# (6)  PHQi9 + baseline variables
nice_lasso_coefs |> 
  filter(Varset %in% c(6, 11, 15, 16, 17)) |> 
  group_by(Varset, Variable) |> 
  summarize(across(Estimate, .fns = list(mn = ~ mean(.x), sd = ~ sd(.x), 
                                         min = ~ min(.x), max = ~ max(.x),
                                         pct25 = ~ quantile(.x, .25),
                                         pct50 = ~ median(.x),
                                         pct75 = ~ quantile(.x, .75)),
                   .names = "{.fn}_{.col}")) |> 
  readr::write_csv(file = paste0(results_dir, "summ_lasso_coefs.csv"))

phq9_lasso_plot <- nice_lasso_coefs |> 
  filter(Varset %in% c(6, 16, 17), grepl("item9", Variable)) |> 
  ggplot(aes(x = Timepoint, y = Estimate, color = nice_varset, shape = nice_varset)) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(color = "Variable group", shape = "Variable group") 
ggsave(filename = paste0(results_dir, "plot_lasso_phq9.png"), width = 10, height = 5,
       plot = phq9_lasso_plot)
