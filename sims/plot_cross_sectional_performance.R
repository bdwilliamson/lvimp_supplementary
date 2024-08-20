# load simulation results, create plots

# load required functions and packages -----------------------------------------
library("tidyr")
library("dplyr")
library("data.table")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("here")

source(here::here("sims", "utils.R"))

if (!dir.exists(here::here("..", "plots", "sims"))) {
  dir.create(here::here("..", "plots", "sims"), recursive = TRUE)
}
plots_dir <- paste0(here::here("..", "plots", "sims"), "/")
results_dir <- here::here("..", "results", "sims")
results_dir <- "H:/Papers/longitudinal_vim/results/sims/"

# read in results --------------------------------------------------------------
# these results come from compile_cross_sectional_performance.R
all_output <- readRDS(paste0(results_dir, "all_output.rds"))

# create plots -----------------------------------------------------------------
# since the simulations for the true values are 2500 replications, can go to 3rd decimal
round_digits <- 3
output_tib <- as_tibble(all_output) %>% 
  mutate(truth = round(truth, round_digits)) %>% 
  mutate(bias = (est - truth), cover_init = cil <= truth & ciu >= truth,
         width_init = ciu - cil, reject_init = p_value < 0.05,
         estimator = factor(algo, labels = c("GLM", "LASSO", "RF", "SL", "XGB"),
                            levels = c("glm", "glmnet", "rf", "SL", "xgb")),
         varset_fct = factor(varset))
summary_tib <- output_tib %>% 
  group_by(n, estimator, outcome_type, corr_between, corr_within, varset_fct,
           designation) %>% 
  summarize(mn_est = mean(est), mdn_est = median(est),
            truth = mean(truth),
            bias = mean(bias), mdn_bias = median(bias),
            ese = sd(est),
            cover = mean(cover_init), reject = mean(reject_init, na.rm = TRUE),
            width = mean(width_init), .groups = "drop")

all_varsets <- unique(output_tib$varset)
all_designations <- unique(output_tib$designation)
all_corrs <- unique(output_tib$corr_within)

# plots of bias, coverage, power
dodge_width <- 0.875
legend_text_size <- 10
axis_text_size <- 10
title_text_size <- 12
fig_width <- 9
fig_height <- 6
point_size <- 1.5
label_x <- 0.025
bias_ylim <- c(-0.45, 0.45)
reject_ylim <- c(0, 1)
cover_ylim <- c(0.5, 1)
width_ylim <- c(0, 0.5)
for (k in seq_len(length(all_corrs))) {
  this_corr <- all_corrs[k]
  for (i in seq_len(length(all_varsets))) {
    this_varset <- all_varsets[i]
    for (j in seq_len(length(all_designations))) {
      this_desig <- all_designations[j]
      if (grepl("vim", this_desig) & (grepl("baseline", this_varset) | grepl("all", this_varset))) {
        # do nothing
      } else {
        # plot name suffix
        plot_name_suffix <- nice_vim_designation(vim_descr = this_varset,
                                                 designation = this_desig)
        # subset
        this_output_tib <- output_tib %>% 
          filter(varset == this_varset, designation == this_desig,
                 corr_within == this_corr)
        this_summ_tib <- summary_tib %>% 
          filter(varset_fct == this_varset, designation == this_desig,
                 corr_within == this_corr)
        shapes <- c(16, 17, 15, 3, 7)
        if (this_corr == 0.5) {
          this_output_tib <- this_output_tib %>% 
            filter(estimator == "GLM")
          this_summ_tib <- this_summ_tib %>% 
            filter(estimator == "GLM")
          shapes <- 16
        }
        # plot bias
        bias_plot <- this_output_tib %>% 
          ggplot(aes(x = estimator, y = bias, color = estimator)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          geom_boxplot(position = position_dodge(dodge_width)) +
          scale_color_viridis_d(begin = 0, end = 0.75) +
          ylab(expression(paste("empirical ", bias[n], sep = ""))) +
          ylim(bias_ylim) +
          xlab("n") +
          labs(shape = "Estimator", color = "Estimator") +
          ggtitle(paste0("BIAS")) +
          facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
          theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
                panel.spacing = unit(0, "cm"),
                strip.background = element_blank(), strip.placement = "outside",
                panel.grid.minor.x = element_line(color = "grey85"),
                panel.grid.major.y = element_line(color = "grey85")) +
          geom_vline(aes(xintercept = 0.4), color = "grey85")
        # plot coverage
        cover_plot <- this_summ_tib %>% 
          ggplot(aes(x = estimator, y = cover, shape = estimator, color = estimator)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
          geom_point(position = position_dodge(dodge_width), size = point_size) +
          scale_color_viridis_d(begin = 0, end = 0.75) +
          ylab("Empirical coverage") +
          ylim(cover_ylim) +
          xlab("n") +
          labs(shape = "Estimator", color = "Estimator") +
          ggtitle(paste0("COVERAGE")) +
          facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
          theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
                panel.spacing = unit(0, "cm"),
                strip.background = element_blank(), strip.placement = "outside",
                panel.grid.minor.x = element_line(color = "grey85"),
                panel.grid.major.y = element_line(color = "grey85")) +
          geom_vline(aes(xintercept = 0.4), color = "grey85")
        # plot CI width
        width_plot <- this_summ_tib %>% 
          ggplot(aes(x = estimator, y = width, shape = estimator, color = estimator)) +
          geom_point(position = position_dodge(dodge_width), size = point_size) +
          scale_color_viridis_d(begin = 0, end = 0.75) +
          ylab("Confidence interval width") +
          ylim(width_ylim) +
          xlab("n") +
          labs(shape = "Estimator", color = "Estimator") +
          ggtitle(paste0("WIDTH")) +
          facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
          theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
                panel.spacing = unit(0, "cm"),
                strip.background = element_blank(), strip.placement = "outside",
                panel.grid.minor.x = element_line(color = "grey85"),
                panel.grid.major.y = element_line(color = "grey85")) +
          geom_vline(aes(xintercept = 0.4), color = "grey85")
        # plot power/type I error
        power_plot <- this_summ_tib %>% 
          ggplot(aes(x = estimator, y = reject, shape = estimator, color = estimator)) +
          geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
          geom_point(position = position_dodge(dodge_width), size = point_size) +
          scale_color_viridis_d(begin = 0, end = 0.75) +
          ylab("Proportion of tests rejected") +
          ylim(reject_ylim) +
          xlab("n") +
          labs(shape = "Estimator", color = "Estimator") +
          ggtitle(paste0("REJECTION PROPORTION")) +
          facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
          theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
                panel.spacing = unit(0, "cm"),
                strip.background = element_blank(), strip.placement = "outside",
                panel.grid.minor.x = element_line(color = "grey85"),
                panel.grid.major.y = element_line(color = "grey85")) +
          geom_vline(aes(xintercept = 0.4), color = "grey85")
        # combine the plots
        summ_for_legend <- mutate(this_summ_tib, est = estimator)
        common_legend <- get_legend(
          bias_plot +
            guides(color = guide_legend(nrow = 2),
                   shape = guide_legend(nrow = 2, override.aes = list(shape = rep(NA, length(unique(summ_for_legend$est)))))) +
            geom_point(aes(x = estimator, y = cover, shape = est, alpha = est), data = summ_for_legend) +
            scale_alpha_manual(name = NULL, values = rep(1, length(unique(summ_for_legend$est))),
                               breaks = unique(summ_for_legend$est),
                               guide = guide_legend(nrow = 2,
                                                    override.aes = list(
                                                      color = "black",
                                                      shape = shapes
                                                    ))) +
            theme(legend.direction = "horizontal",
                  legend.position = "bottom",
                  legend.title = element_text(size = legend_text_size),
                  legend.text = element_text(size = legend_text_size),
                  legend.spacing.x = unit(0.5, "cm"))
        ) 
        if (grepl("baseline", this_varset)) {
          combined_plot <- plot_grid(
            bias_plot + theme(legend.position = "none",
                              title = element_text(size = title_text_size),
                              axis.title = element_text(size = axis_text_size),
                              axis.text = element_text(size = axis_text_size),
                              plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            cover_plot + theme(legend.position = "none",
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            width_plot + theme(legend.position = "none",
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            labels = "AUTO", label_x = label_x, label_size = title_text_size
          )
        } else {
          combined_plot <- plot_grid(
            bias_plot + theme(legend.position = "none",
                              title = element_text(size = title_text_size),
                              axis.title = element_text(size = axis_text_size),
                              axis.text = element_text(size = axis_text_size),
                              plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            cover_plot + theme(legend.position = "none",
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            power_plot + theme(legend.position = "none",
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            width_plot + theme(legend.position = "none",
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
            labels = "AUTO", label_x = label_x, label_size = title_text_size
          )
        }
        final_plot <- plot_grid(combined_plot, common_legend, ncol = 1, nrow = 2, 
                                rel_heights = c(1, .1))
        # save off the plot
        ggsave(filename = paste0(plots_dir, "/individual_varsets/varset_", this_varset, "_", this_desig, "_corr_", this_corr, ".png"),
               final_plot, width = fig_width, height = fig_height, units = "in",
               dpi = 300) 
      }
    }
  }
}

# summary plots and tables for main manuscript ---------------------------------
varsets_of_interest <- expand.grid(type = c("addi_", "loco_"), varset = c(1, 2, 3, 8))
varsets_of_interest$raw_varset <- paste0(varsets_of_interest$type, varsets_of_interest$varset)
varsets_of_interest$nice_varset <- gsub("loco_", "Leave-out: ", gsub("addi_", "Add-in: ", varsets_of_interest$raw_varset))

# true values
true_values <- output_tib |> 
  group_by(p, outcome_type, corr_between, corr_within, varset, designation) |> 
  slice(1) |> 
  select(varset, designation, truth)

nice_true_values <- true_values |> 
  filter(varset %in% varsets_of_interest$raw_varset, grepl("vim", designation) & !grepl("autc", designation)) |> 
  mutate(varset_fct = factor(case_when(
    varset == "addi_1" ~ "Add-in: 1",
    varset == "addi_2" ~ "Add-in: 2",
    varset == "addi_3" ~ "Add-in: 3",
    varset == "addi_8" ~ "Add-in: 8",
    varset == "loco_1" ~ "Leave-out: 1",
    varset == "loco_2" ~ "Leave-out: 2",
    varset == "loco_3" ~ "Leave-out: 3",
    varset == "loco_8" ~ "Leave-out: 8"
  ), levels = varsets_of_interest$nice_varset, ordered = TRUE)) |> 
  ungroup() |> 
  mutate(designation = gsub("vim-", "", designation),
         designation = ifelse(!is.na(as.numeric(designation)), paste0("Timepoint ", designation), designation)) |> 
  arrange(designation, varset_fct) |> 
  select(designation, varset_fct, truth, corr_within) |> 
  mutate(vim_type = ifelse(grepl("Add-in", varset_fct), "Add-in", "Leave-out"),
         variable = gsub("Add-in: ", "", gsub("Leave-out: ", "", varset_fct)))

wide_true_values <- nice_true_values |> 
  select(-vim_type, -variable) |> 
  pivot_wider(names_from = designation, values_from = truth) |> 
  rename(`VIM type: variable` = varset_fct, `Mean` = average, 
         `Trend intercept` = `trend-intercept`, `Trend slope` = `trend-slope`) |> 
  mutate(across(where(is.numeric), ~ as.character(round(.x, digits = 3))))

# a plot with true variable importance values (and summaries)
summary_only_table <- nice_true_values |> 
  filter(!grepl("Timepoint", designation)) |> 
  mutate(designation = case_when(
    designation == "average" ~ "Mean",
    designation == "trend-intercept" ~ "Trend: intercept",
    designation == "trend-slope" ~ "Trend: slope"
  ))
font_size <- 9

for (i in 1:length(all_corrs)) {
  this_corr <- all_corrs[i]
  true_vims_plot <- ggplot(data = nice_true_values |> 
                               filter(grepl("Timepoint", designation), corr_within == this_corr) |> 
                               mutate(designation = gsub("Timepoint ", "", designation)), 
                             aes(x = designation, y = truth, shape = vim_type)) +
    geom_point(position = position_dodge(width = 0.2)) + 
    scale_color_viridis_d(begin = 0, end = 0.75) +
    labs(shape = "VIM type", x = "Timepoint", y = "True VIM value") +
    facet_grid(cols = vars(variable), labeller = label_both) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  true_summaries_plot <- ggplot(data = nice_true_values |>
                                    filter(!grepl("Timepoint", designation), corr_within == this_corr) |>
                                    mutate(designation = case_when(
                                      designation == "average" ~ "Mean",
                                      designation == "trend-intercept" ~ "Trend: intercept",
                                      designation == "trend-slope" ~ "Trend: slope"
                                    ),
                                    xval = 1,
                                    yval = case_when(
                                      designation == "Mean" ~ 0.01,
                                      designation == "Trend: intercept" ~ 0.05,
                                      designation == "Trend: slope" ~ 0.03
                                    ))) +
    geom_text(aes(x = xval, y = yval, label = paste0(designation, " = ", round(truth, 3))),
              hjust = 0, size = 4) +
    xlim(c(1, 4)) +
    ylim(c(0, 0.06)) +
    facet_grid(cols = vars(variable), rows = vars(vim_type)) +
    theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0.75), units = "in")) 
  true_vim_plot <- cowplot::plot_grid(
    true_vims_plot, true_summaries_plot, 
    rel_heights = c(1, 0.75), nrow = 2, ncol = 1
  )  
  ggsave(filename = paste0(plots_dir, "sim_true_vims_", this_corr, ".png"), true_vim_plot,
         width = 8.5, height = 4, units = "in")
  # a table with true variable importance values
  corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
  knitr::kable(wide_true_values %>% filter(corr_within == this_corr) %>% select(-corr_within), 
               digits = 3, format = "latex", booktabs = TRUE,
               caption = paste0("True variable importance values at each time point and summarized over the time series.",
                                corr_txt,
                                " \\label{tab:sim_true_vims_", this_corr, "}"),
               linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size) |> 
    kableExtra::save_kable(file = paste0(plots_dir, "sim_true_vims_", this_corr, ".tex"))
  
}

  
# subset
plot_tib <- summary_tib %>% 
  filter(varset_fct %in% varsets_of_interest$raw_varset, designation %in% paste0("vim-", 1:4)) %>% 
  mutate(timepoint = gsub("vim-", "", designation)) %>% 
  mutate(varset_fct = factor(case_when(
    varset_fct == "addi_1" ~ "Add-in: 1",
    varset_fct == "addi_2" ~ "Add-in: 2",
    varset_fct == "addi_3" ~ "Add-in: 3",
    varset_fct == "addi_8" ~ "Add-in: 8",
    varset_fct == "loco_1" ~ "Leave-out: 1",
    varset_fct == "loco_2" ~ "Leave-out: 2",
    varset_fct == "loco_3" ~ "Leave-out: 3",
    varset_fct == "loco_8" ~ "Leave-out: 8"
  ), levels = varsets_of_interest$nice_varset, ordered = TRUE), 
  n_fct = factor(paste0("n: ", n), levels = paste0("n: ", c(100, 250, 500, 1000, 5000, 10000)))) |> 
  mutate(vim_type = factor(ifelse(grepl("add", varset_fct, ignore.case = TRUE), "Add-in", "Leave-out")),
         variable = factor(
           paste0("variable: ", gsub("Add-in: ", "", gsub("Leave-out: ", "", varset_fct)))
         ))
summ_tib <- summary_tib %>% 
  filter(varset_fct %in% varsets_of_interest$raw_varset, designation %in% c("vim-average", "vim-trend-slope")) %>% 
  mutate(varset_fct = factor(case_when(
    varset_fct == "addi_1" ~ "Add-in: 1",
    varset_fct == "addi_2" ~ "Add-in: 2",
    varset_fct == "addi_3" ~ "Add-in: 3",
    varset_fct == "addi_8" ~ "Add-in: 8",
    varset_fct == "loco_1" ~ "Leave-out: 1",
    varset_fct == "loco_2" ~ "Leave-out: 2",
    varset_fct == "loco_3" ~ "Leave-out: 3",
    varset_fct == "loco_8" ~ "Leave-out: 8"
  ), levels = varsets_of_interest$nice_varset, ordered = TRUE)) |> 
  mutate(designation = gsub("vim-", "", designation)) |> 
  arrange(estimator, designation, varset_fct) |> 
  select(designation, varset_fct, estimator, everything(), -bias, -mdn_est, -mdn_bias) |> 
  rename(Estimator = estimator, `VIM type: variable` = varset_fct,
         `$n$` = n,
         `Summary` = designation, `Mean est.` = mn_est, `True value` = truth, 
         `Empirical SE` = ese, Coverage = cover, `Rejection prop.` = reject, `CI width` = width) |> 
  mutate(across(.cols = where(is.numeric), ~ round(.x, digits = 3))) |> 
  mutate(across(.cols = where(is.numeric), ~ as.character(.x)))

# plot of add-in  & leave-out variable importance over time for variables 1, 2, 3, 8
# make it look like the figure for true VIMs
for (i in 1:length(all_corrs)) {
  this_corr <- all_corrs[i]
  this_plot_tib <- plot_tib %>% 
    filter(corr_within == this_corr)
  if (this_corr == 0.5) {
    this_plot_tib <- this_plot_tib %>% 
      filter(estimator == "GLM")
  }
  summary_plot_sl_glm <- this_plot_tib %>% 
    filter(estimator %in% c("GLM", "SL")) |> 
    ggplot(aes(x = timepoint, y = mn_est, color = estimator, shape = vim_type)) +
    geom_point(position = position_dodge(dodge_width)) +
    scale_color_viridis_d(begin = 0, end = 0.75) +
    labs(x = "Timepoint", y = "Mean value of estimated VIM",
         color = "Estimator", shape = "VIM type") +
    ylim(c(-0.01, 0.35)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_grid(rows = vars(variable), cols = vars(n_fct))
  summary_plot_minus_sl_glm <- this_plot_tib %>% 
    filter(!(estimator %in% c("GLM", "SL"))) |> 
    ggplot(aes(x = timepoint, y = mn_est, color = estimator, shape = vim_type)) +
    geom_point(position = position_dodge(dodge_width)) +
    scale_color_viridis_d(begin = 0, end = 0.75) +
    labs(x = "Timepoint", y = "Mean value of estimated VIM",
         color = "Estimator", shape = "VIM type") +
    ylim(c(-0.01, 0.35)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_grid(rows = vars(variable), cols = vars(n_fct))
  ggsave(filename = paste0(plots_dir, "vims_of_interest_over_time_", this_corr, ".png"),
         summary_plot_sl_glm, width = fig_width, height = fig_height, units = "in",
         dpi = 300)
  if (this_corr == 0) {
    ggsave(filename = paste0(plots_dir, "vims_of_interest_over_time_supp_", this_corr, ".png"),
           summary_plot_minus_sl_glm, width = fig_width, height = fig_height, units = "in",
           dpi = 300)  
  }
}

# table of variable importance for linear trend, average (sample size in meta-rows)
# columns = variable, summary, importance type, bias, variance, coverage, power
# pick one (large) n
font_size <- 9
for (i in 1:length(all_corrs)) {
  this_corr <- all_corrs[i]
  corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
  this_summ_tib <- summ_tib %>% 
    filter(corr_within == this_corr)
  all_est_txt <- paste0("all estimators: logistic regression (GLM),",
                        " lasso, random forests (RF), Super Learner (SL), and boosted trees (XGB)."
  )
  if (this_corr == 0.5) {
    this_summ_tib <- this_summ_tib %>% 
      filter(Estimator == "GLM")
    all_est_txt <- paste0("logistic regression (GLM).")
  }
  summ_tab_bign <- this_summ_tib |> 
    filter(`$n$` == 5000) |> 
    select(-outcome_type, -corr_between, -corr_within, -`$n$`)
  knitr::kable(summ_tab_bign |> 
                 filter(grepl("average", Summary)), format = "latex", booktabs = TRUE,
               caption = paste0("Performance of estimators of the average VIM over the time series,",
                                " summarized over 1000 Monte-Carlo replications. Performance",
                                " at sample size $n = 5000$ shown for ", all_est_txt,
                                " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                " variable importance are displayed for each variable.", corr_txt,
                                " \\label{tab:sim_bign_avg_", this_corr, "}"),
               linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size) |> 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_avg_", this_corr, ".tex"))
  
  knitr::kable(summ_tab_bign |> 
                 filter(grepl("slope", Summary)), format = "latex", booktabs = TRUE,
               caption = paste0("Performance of estimators of the slope of the linear VIM trend over the time series,",
                                " summarized over 1000 Monte-Carlo replications. Performance",
                                " at sample size $n = 5000$ shown for ", all_est_txt, 
                                " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                " variable importance are displayed for each variable.", corr_txt,
                                " \\label{tab:sim_bign_slope_", this_corr, "}"),
               linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size) |> 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_slope_", this_corr, ".tex"))
  
  # pick one estimator (SL)
  summ_tab_sl <- summ_tib |> 
    filter(Estimator == "SL") |> 
    select(-outcome_type, -corr_between, -corr_within, -Estimator)
  knitr::kable(summ_tab_sl, format = "latex", booktabs = TRUE, escape = FALSE,
               caption = paste0("Performance of estimators of the average VIM and ",
                                " the slope of the linear VIM trend over the time series,",
                                " summarized over 1000 Monte-Carlo replications. Performance",
                                " shown for each sample size $n \\in \\{100, 250, 1000, 5000, 10000\\}$,",
                                " using the Super Learner.",
                                " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                " variable importance are displayed for each variable.", corr_txt,
                                " \\label{tab:sim_sl_", this_corr, "}"),
               linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size) |> 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_sl_", this_corr, ".tex"))
  
  # separate average from slope
  summ_tib |> 
    filter(grepl("average", `Summary`)) |> 
    select(-outcome_type, -corr_between, -corr_within, -`Summary`) |> 
    knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                 caption = paste0("Performance of estimators of the average VIM over the time series,",
                                  " summarized over 1000 Monte-Carlo replications. Performance",
                                  " shown for each sample size $n \\in \\{100, 250, 1000, 5000, 10000\\}$" ,
                                  " and ", all_est_txt,
                                  " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                  " variable importance are displayed for each variable.", corr_txt,
                                  " \\label{tab:sim_all_avg_", this_corr, "}"),
                 linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size, 
                              latex_options = c("repeat_header"),
                              repeat_header_text = paste0("Performance of estimators of the average VIM over the time series,",
                                                          " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                              repeat_header_method = "replace") |> 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_avg_", this_corr, ".tex"))
  summ_tib |> 
    filter(grepl("slope", `Summary`)) |> 
    select(-outcome_type, -corr_between, -corr_within, -`Summary`) |> 
    knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                 caption = paste0("Performance of estimators of the linear VIM trend over the time series,",
                                  " summarized over 1000 Monte-Carlo replications. Performance",
                                  " shown for each sample size $n \\in \\{100, 250, 1000, 5000, 10000\\}$" ,
                                  " and ", all_est_txt, 
                                  " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                  " variable importance are displayed for each variable.", corr_txt,
                                  " \\label{tab:sim_all_slope_", this_corr, "}"),
                 linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size, 
                              latex_options = c("repeat_header"),
                              repeat_header_text = paste0("Performance of estimators of the linear VIM trend over the time series,",
                                                          " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                              repeat_header_method = "replace") |> 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_slope_", this_corr, ".tex"))
}


# also print out separated by X to make the table less long
for (j in 1:length(all_corrs)) {
  this_corr <- all_corrs[j]
  corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
  for (i in c(1, 2, 3, 8)) {
    this_summ_tib <- summ_tib |> 
      filter(grepl(i, summ_tib$`VIM type: variable`), corr_within == this_corr)
    all_est_txt <- paste0("all estimators: logistic regression (GLM),",
                          " lasso, random forests (RF), Super Learner (SL), and boosted trees (XGB)."
    )
    if (this_corr == 0.5) {
      this_summ_tib <- this_summ_tib %>% 
        filter(Estimator == "GLM")
      all_est_txt <- "logistic regression (GLM)."
    }
    this_summ_tib |> 
      filter(grepl("average", `Summary`)) |> 
      select(-outcome_type, -corr_between, -corr_within, -`Summary`) |> 
      knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                   caption = paste0("Performance of estimators of the average VIM for variable ", i, " over the time series,",
                                    " summarized over 1000 Monte-Carlo replications. Performance",
                                    " shown for each sample size $n \\in \\{100, 250, 1000, 5000, 10000\\}$" ,
                                    " and ", all_est_txt, 
                                    " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                    " variable importance are displayed.", corr_txt,
                                    " \\label{tab:sim_all_avg_", i, "_", this_corr, "}"),
                   linesep = "") |> 
      kableExtra::kable_styling(font_size = font_size, 
                                latex_options = c("repeat_header"),
                                repeat_header_text = paste0("Performance of estimators of the average VIM for variable ", i, " over the time series,",
                                                            " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                repeat_header_method = "replace") |> 
      kableExtra::save_kable(file = paste0(plots_dir, "summary_table_avg_", i, "_", this_corr, ".tex"))
    this_summ_tib |> 
      filter(grepl("slope", `Summary`)) |> 
      select(-outcome_type, -corr_between, -corr_within, -`Summary`) |> 
      knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                   caption = paste0("Performance of estimators of the linear VIM trend for variable ", i, " over the time series,",
                                    " summarized over 1000 Monte-Carlo replications. Performance",
                                    " shown for each sample size $n \\in \\{100, 250, 1000, 5000, 10000\\}$" ,
                                    " and ", all_est_txt, 
                                    " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                    " variable importance are displayed.", corr_txt,
                                    " \\label{tab:sim_all_slope_", i, "_", this_corr, "}"),
                   linesep = "") |> 
      kableExtra::kable_styling(font_size = font_size, 
                                latex_options = c("repeat_header"),
                                repeat_header_text = paste0("Performance of estimators of the linear VIM trend for variable ", i, " over the time series,",
                                                            " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                repeat_header_method = "replace") |> 
      kableExtra::save_kable(file = paste0(plots_dir, "summary_table_slope_", i, "_", this_corr, ".tex"))
  }
  
}

summ_tab_bign <- summ_tib |> 
  # filter((`$n$` == 5000 & corr_within == 0) | (`$n$` == 1000 & corr_within == 0.5)) |> 
  filter(`$n$` == 5000) %>% 
  select(-outcome_type, -corr_between, -`$n$`)


# summarize only for SL, GLM at n = 5000
for (i in 1:length(all_corrs)) {
  this_corr <- all_corrs[i]
  corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
  all_est_txt <- paste0("logistic regression (GLM) and Super Learner (SL).")
  n_txt <- "$n = 5000$"
  this_summ_tab_bign <- summ_tab_bign %>% 
    filter(corr_within == this_corr) %>% 
    select(-corr_within)
  
  if (this_corr == 0.5) {
    this_summ_tab_bign <- this_summ_tab_bign %>% 
      filter(Estimator == "GLM")
    all_est_txt <- "logistic regression (GLM)."
    # n_txt <- "$n = 1000$"
  } 
  big_n_summ <- this_summ_tab_bign |> 
    filter(Estimator == "GLM" | Estimator == "SL") |> 
    select(-Summary, -Estimator) |> 
    knitr::kable(format = "latex", booktabs = TRUE,
                 caption = paste0("Performance of estimators of the average VIM and slope of the linear trend in VIM over the time series,",
                                  " summarized over 1000 Monte-Carlo replications. Performance",
                                  " at sample size ", n_txt, " shown for ", all_est_txt, 
                                  " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                  " variable importance are displayed for each variable. The true importance of variable 1",
                                  " is nearly constant over time; variable 2 has increasing importance over time; variable 3",
                                  " has decreasing importance over time; and variable 8 has zero importance at all time points.", corr_txt,
                                  " \\label{tab:sim_bign_all_", this_corr, "}"),
                 linesep = "") |> 
    kableExtra::kable_styling(font_size = font_size, 
                              latex_options = c("scale_down")) |> 
    kableExtra::pack_rows("Average (GLM)", 1, 8) |> 
    kableExtra::pack_rows("Trend - slope (GLM)", 9, 16)
  if (this_corr == 0) {
    big_n_summ <- big_n_summ %>% 
      kableExtra::pack_rows("Average (SL)", 17, 24) |> 
      kableExtra::pack_rows("Trend - slope (SL)", 25, 32)
  }
  big_n_summ %>% 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_all_", this_corr, ".tex"))
}
