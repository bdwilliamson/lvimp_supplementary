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

# read in results --------------------------------------------------------------
# these results come from compile_cross_sectional_performance.R
all_output <- readRDS(paste0(results_dir, "all_output.rds"))

# create plots -----------------------------------------------------------------
# since the simulations for the true values are 2500 replications, can go to 3rd decimal
round_digits <- 3
output_tib <- all_output_tib %>% 
  mutate(truth = round(truth, round_digits)) %>% 
  mutate(bias = (est - truth), cover_init = cil <= truth & ciu >= truth,
         width_init = ciu - cil, reject_init = p_value < 0.05,
         estimator = factor(algo, labels = c("GLM", "LASSO", "RF", "SL", "XGB"),
                            levels = c("glm", "glmnet", "rf", "SL", "xgb")),
         varset_fct = factor(varset))
eses <- output_tib %>% 
  group_by(measure, n, estimator, outcome_type, corr_between, corr_within, varset_fct,
           designation, dgm) %>% 
  summarize(ese = sd(est, na.rm = TRUE))
summary_tib <- output_tib %>% 
  left_join(eses, by = c("measure", "n", "estimator", "outcome_type", "corr_between", 
                         "corr_within", "varset_fct",
                         "designation", "dgm")) %>% 
  mutate(ese_cover_init = est - qnorm(0.975) * ese <= truth & truth <= est + qnorm(0.975) * ese) %>% 
  group_by(measure, n, estimator, outcome_type, corr_between, corr_within, varset_fct,
           designation, dgm) %>% 
  summarize(mn_est = mean(est, na.rm = TRUE), mdn_est = median(est, na.rm = TRUE),
            truth = mean(truth, na.rm = TRUE),
            bias = mean(bias, na.rm = TRUE), mdn_bias = median(bias, na.rm = TRUE),
            ese = mean(ese, na.rm = TRUE),
            cover = mean(cover_init, na.rm = TRUE), 
            ese_cover = mean(ese_cover_init, na.rm = TRUE),
            reject = mean(reject_init, na.rm = TRUE),
            width = mean(width_init, na.rm = TRUE), .groups = "drop") %>% 
  mutate(nice_measure = case_when(
    measure == "auc" | measure == "ppv" ~ toupper(measure),
    measure == "sensitivity" ~ "Sensitivity"
  )) %>% 
  select(-ese_cover)

all_varsets <- unique(output_tib$varset)
all_designations <- unique(output_tib$designation)
all_corrs <- unique(output_tib$corr_within)
all_measures <- unique(output_tib$measure)
all_dgms <- unique(output_tib$dgm)
all_nice_measures <- case_when(
  all_measures == "auc" | all_measures == "ppv" ~ toupper(all_measures),
  all_measures == "sensitivity" ~ "Sensitivity"
)

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

for (d in seq_len(length(all_dgms))) {
  this_dgm <- all_dgms[d]
  for (l in seq_len(length(all_measures))) {
    this_measure <- all_measures[l]
    for (k in seq_len(length(all_corrs))) {
      this_corr <- all_corrs[k]
      for (i in seq_len(length(all_varsets))) {
        this_varset <- all_varsets[i]
        for (j in seq_len(length(all_designations))) {
          this_desig <- all_designations[j]
          this_output_tib <- output_tib %>% 
            filter(varset == this_varset, designation == this_desig,
                   corr_within == this_corr, measure == this_measure,
                   dgm == this_dgm)
          max_bias <- max(abs(this_output_tib$bias))
          bias_ylim <- c(-max_bias, max_bias)
          if (grepl("vim", this_desig) & (grepl("baseline", this_varset) | grepl("all", this_varset))) {
            # do nothing
          } else if (nrow(this_output_tib) == 0) {
            # also do nothing
          } else {
            # plot name suffix
            plot_name_suffix <- nice_vim_designation(vim_descr = this_varset,
                                                     designation = this_desig)
            # subset
            this_summ_tib <- summary_tib %>% 
              filter(varset_fct == this_varset, designation == this_desig,
                     corr_within == this_corr, measure == this_measure)
            shapes <- c(16, 17, 15, 3, 7)
            if (this_corr == 0.5) {
              this_output_tib <- this_output_tib %>%
                filter(estimator == "SL")
              this_summ_tib <- this_summ_tib %>%
                filter(estimator == "SL")
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
            ggsave(filename = paste0(plots_dir, "/individal_varsets/varset_", this_varset, "_", this_desig, "_corr_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".png"),
                   final_plot, width = fig_width, height = fig_height, units = "in",
                   dpi = 300) 
          }
        }
      }
    }
  }  
}


# summary plots and tables for main manuscript ---------------------------------
varsets_of_interest <- expand.grid(type = c("addi_", "loco_"), varset = c(1, 2, 3, 8))
varsets_of_interest$raw_varset <- paste0(varsets_of_interest$type, varsets_of_interest$varset)
varsets_of_interest$nice_varset <- gsub("loco_", "Leave-out: ", gsub("addi_", "Add-in: ", varsets_of_interest$raw_varset))

# true values
true_values <- output_tib %>% 
  group_by(p, outcome_type, measure, corr_between, corr_within, varset, designation, dgm) %>% 
  slice(1) %>% 
  select(varset, designation, dgm, truth)

nice_true_values <- true_values %>% 
  filter(varset %in% varsets_of_interest$raw_varset, grepl("vim", designation) & !grepl("autc", designation)) %>% 
  mutate(varset_fct = factor(case_when(
    varset == "addi_1" ~ "Add-in: 1",
    varset == "addi_2" ~ "Add-in: 2",
    varset == "addi_3" ~ "Add-in: 3",
    varset == "addi_8" ~ "Add-in: 8",
    varset == "loco_1" ~ "Leave-out: 1",
    varset == "loco_2" ~ "Leave-out: 2",
    varset == "loco_3" ~ "Leave-out: 3",
    varset == "loco_8" ~ "Leave-out: 8"
  ), levels = varsets_of_interest$nice_varset, ordered = TRUE)) %>% 
  ungroup() %>% 
  mutate(designation = gsub("vim-", "", designation),
         designation = ifelse(!is.na(as.numeric(designation)), paste0("Timepoint ", designation), designation)) %>% 
  arrange(designation, varset_fct) %>% 
  select(measure, designation, dgm, varset_fct, truth, corr_within) %>% 
  mutate(vim_type = ifelse(grepl("Add-in", varset_fct), "Add-in", "Leave-out"),
         variable = gsub("Add-in: ", "", gsub("Leave-out: ", "", varset_fct))) %>% 
  mutate(nice_measure = case_when(
    measure == "auc" | measure == "ppv" ~ toupper(measure),
    measure == "sensitivity" ~ "Sensitivity"
  ))

wide_true_values <- nice_true_values %>% 
  select(-vim_type, -variable) %>% 
  pivot_wider(names_from = designation, values_from = truth) %>% 
  rename(`VIM type: variable` = varset_fct, `Mean` = average, 
         `Trend intercept` = `trend-intercept`, `Trend slope` = `trend-slope`) %>% 
  mutate(across(where(is.numeric), ~ as.character(round(.x, digits = 3))))

nice_true_values_pred <- true_values %>% 
  filter(varset %in% c(varsets_of_interest$raw_varset, "baseline", "all"), grepl("predictiveness", designation) & !grepl("autc", designation)) %>% 
  mutate(varset_fct = factor(case_when(
    varset == "addi_1" ~ "Add-in: 1",
    varset == "addi_2" ~ "Add-in: 2",
    varset == "addi_3" ~ "Add-in: 3",
    varset == "addi_8" ~ "Add-in: 8",
    varset == "loco_1" ~ "Leave-out: 1",
    varset == "loco_2" ~ "Leave-out: 2",
    varset == "loco_3" ~ "Leave-out: 3",
    varset == "loco_8" ~ "Leave-out: 8",
    varset == "baseline" ~ "Baseline",
    varset == "all" ~ "All"
  ), levels = c(varsets_of_interest$nice_varset, "Baseline", "All"), ordered = TRUE)) %>% 
  ungroup() %>% 
  mutate(designation = gsub("predictiveness-", "", designation),
         designation = ifelse(!is.na(as.numeric(designation)), paste0("Timepoint ", designation), designation)) %>% 
  arrange(designation, varset_fct) %>% 
  select(measure, designation, dgm, varset_fct, truth, corr_within) %>% 
  mutate(vim_type = ifelse(grepl("Add-in", varset_fct), "Add-in", "Leave-out"),
         variable = gsub("Add-in: ", "", gsub("Leave-out: ", "", varset_fct))) %>% 
  mutate(nice_measure = case_when(
    measure == "auc" | measure == "ppv" ~ toupper(measure),
    measure == "sensitivity" ~ "Sensitivity"
  ))

wide_true_values_pred <- nice_true_values_pred %>% 
  select(-vim_type, -variable) %>% 
  pivot_wider(names_from = designation, values_from = truth) %>% 
  rename(`Pred. type: variable` = varset_fct, `Mean` = average, 
         `Trend intercept` = `trend-intercept`, `Trend slope` = `trend-slope`) %>% 
  mutate(across(where(is.numeric), ~ as.character(round(.x, digits = 3))))


# a plot with true variable importance values (and summaries)
summary_only_table <- nice_true_values %>% 
  filter(!grepl("Timepoint", designation)) %>% 
  mutate(designation = case_when(
    designation == "average" ~ "Mean",
    designation == "trend-intercept" ~ "Trend: intercept",
    designation == "trend-slope" ~ "Trend: slope"
  ))
font_size <- 9

for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  for (i in 1:length(all_corrs)) {
    this_corr <- all_corrs[i]
    for (j in 1:length(all_measures)) {
      this_measure <- all_measures[j]
      this_dataset <- nice_true_values %>% 
        filter(corr_within == this_corr, measure == this_measure, dgm == this_dgm)
      if (nrow(this_dataset) == 0) {
        # do nothing
      } else {
        true_vims_plot <- ggplot(data = this_dataset %>% 
                                   filter(grepl("Timepoint", designation)) %>% 
                                   mutate(designation = gsub("Timepoint ", "", designation)), 
                                 aes(x = designation, y = truth, shape = vim_type)) +
          geom_point(position = position_dodge(width = 0.2)) + 
          scale_color_viridis_d(begin = 0, end = 0.75) +
          labs(shape = "VIM type", x = "Timepoint", y = "True VIM value") +
          facet_grid(cols = vars(variable), labeller = label_both) +
          theme(legend.position = "bottom", legend.direction = "horizontal")
        true_summaries_plot <- ggplot(data = this_dataset %>%
                                        filter(!grepl("Timepoint", designation)) %>%
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
        ggsave(filename = paste0(plots_dir, "sim_true_vims_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".png"), true_vim_plot,
               width = 8.5, height = 4, units = "in")
        # a table with true variable importance values
        corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
        measure_txt <- ifelse(this_measure == "auc", "AUC", 
                              paste0(ifelse(this_measure == "ppv", "PPV", "sensitivity"), " at the 95th percentile of predicted risk"))
        knitr::kable(wide_true_values %>% 
                       filter(corr_within == this_corr, measure == this_measure,
                              dgm == this_dgm) %>% 
                       select(-corr_within, -measure, -dgm, -nice_measure) %>% 
                       mutate(across(where(is.numeric), .fns = ~ as.numeric(sprintf("%.3f", .x)))), 
                     digits = 3, format = "latex", booktabs = TRUE,
                     caption = paste0("True variable importance values (defined using ", measure_txt, ") at each time point and summarized over the time series.",
                                      corr_txt,
                                      " \\label{tab:sim_true_vims_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                     linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size) %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "sim_true_vims_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))   
        # a table with true predictiveness values
        knitr::kable(wide_true_values_pred %>% 
                       filter(corr_within == this_corr, measure == this_measure,
                              dgm == this_dgm) %>% 
                       select(-corr_within, -measure, -dgm, -nice_measure) %>% 
                       mutate(across(where(is.numeric), .fns = ~ as.numeric(sprintf("%.3f", .x)))), 
                     digits = 3, format = "latex", booktabs = TRUE,
                     caption = paste0("True predictiveness values (defined using ", measure_txt, ") at each time point and summarized over the time series.",
                                      corr_txt,
                                      " \\label{tab:sim_true_pred_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                     linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size) %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "sim_true_pred_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))   
      }
    }
  }
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
  n_fct = factor(paste0("n: ", n), levels = paste0("n: ", c(100, 250, 500, 1000, 5000, 10000)))) %>% 
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
  ), levels = varsets_of_interest$nice_varset, ordered = TRUE)) %>% 
  mutate(designation = gsub("vim-", "", designation)) %>% 
  arrange(estimator, designation, varset_fct) %>% 
  select(nice_measure, designation, varset_fct, estimator, everything(), -bias, -mdn_est, -mdn_bias, -measure) %>% 
  rename(Measure = nice_measure, Estimator = estimator, `VIM type: variable` = varset_fct,
         `$n$` = n,
         `Summary` = designation, `Mean est.` = mn_est, `True value` = truth, 
         `Empirical SE` = ese, Coverage = cover, `Rejection prop.` = reject, `CI width` = width) %>% 
  mutate(across(.cols = where(is.numeric), ~ as.numeric(sprintf("%.3f", .x))))
# %>% 
#   mutate(across(.cols = where(is.numeric), ~ as.character(.x)))

# plot of add-in  & leave-out variable importance over time for variables 1, 2, 3, 8
# make it look like the figure for true VIMs
for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  for (i in 1:length(all_corrs)) {
    this_corr <- all_corrs[i]
    for (j in 1:length(all_measures)) {
      this_measure <- all_measures[j]
      this_plot_tib <- plot_tib %>% 
        filter(corr_within == this_corr, measure == this_measure, dgm == this_dgm)
      if (nrow(this_plot_tib) == 0) {
        # do nothing
      } else {
        summary_plot_sl_glm <- this_plot_tib %>% 
          filter(estimator %in% c("GLM", "SL")) %>% 
          ggplot(aes(x = timepoint, y = mn_est, color = estimator, shape = vim_type)) +
          geom_point(position = position_dodge(dodge_width)) +
          scale_color_viridis_d(begin = 0, end = 0.75) +
          labs(x = "Timepoint", y = "Mean value of estimated VIM",
               color = "Estimator", shape = "VIM type") +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
          facet_grid(rows = vars(variable), cols = vars(n_fct))
        ggsave(filename = paste0(plots_dir, "vims_of_interest_over_time_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".png"),
               summary_plot_sl_glm, width = fig_width, height = fig_height, units = "in",
               dpi = 300)
        if ((this_corr == 0 & this_dgm == 1) | (this_corr == 0.5 & this_dgm == 2)) {
          summary_plot_minus_sl_glm <- this_plot_tib %>% 
            filter(!(estimator %in% c("GLM", "SL"))) %>% 
            ggplot(aes(x = timepoint, y = mn_est, color = estimator, shape = vim_type)) +
            geom_point(position = position_dodge(dodge_width)) +
            scale_color_viridis_d(begin = 0, end = 0.75) +
            labs(x = "Timepoint", y = "Mean value of estimated VIM",
                 color = "Estimator", shape = "VIM type") +
            geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
            facet_grid(rows = vars(variable), cols = vars(n_fct))
          ggsave(filename = paste0(plots_dir, "vims_of_interest_over_time_supp_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".png"),
                 summary_plot_minus_sl_glm, width = fig_width, height = fig_height, units = "in",
                 dpi = 300)    
        }
      }
    }
  }
}


# table of variable importance for linear trend, average (sample size in meta-rows)
# columns = variable, summary, importance type, bias, variance, coverage, power
# pick one (large) n
font_size <- 9
for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  for (i in 1:length(all_corrs)) {
    this_corr <- all_corrs[i]
    corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
    for (j in 1:length(all_measures)) {
      this_measure <- all_measures[j]
      this_nice_measure <- all_nice_measures[j]
      measure_txt <- ifelse(this_measure == "auc", "AUC", 
                            paste0(ifelse(this_measure == "ppv", "PPV", "sensitivity"), " at the 95th percentile of predicted risk"))
      this_summ_tib <- summ_tib %>% 
        filter(corr_within == this_corr, Measure == this_nice_measure, dgm == this_dgm) %>% 
        select(-Measure, -dgm)
      if (nrow(this_summ_tib) == 0) {
        # do nothing
      } else {
        all_est_txt <- paste0("all estimators: logistic regression (GLM),",
                              " lasso, random forests (RF), super learner (SL), and boosted trees (XGB)."
        )
        all_measure_txt <- paste0("PPV and sensitivity at the 95th percentile of predicted risk and AUC")
        
        summ_tab_bign <- this_summ_tib %>% 
          filter(`$n$` == 5000) %>% 
          select(-outcome_type, -corr_between, -corr_within, -`$n$`) 
        knitr::kable(summ_tab_bign %>% 
                       filter(grepl("average", Summary)) %>% 
                       select(-Summary), 
                     format = "latex", booktabs = TRUE,
                     caption = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ", over the time series,",
                                      " summarized over 1000 Monte-Carlo replications. Performance",
                                      " at sample size $n = 5000$ shown for ", all_est_txt,
                                      " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                      " variable importance are displayed for each variable.", corr_txt,
                                      " \\label{tab:sim_bign_avg_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                     linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size) %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_avg_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
        
        knitr::kable(summ_tab_bign %>% 
                       filter(grepl("slope", Summary)) %>% 
                       select(-Summary), 
                     format = "latex", booktabs = TRUE,
                     caption = paste0("Performance of estimators of the slope of the linear VIM trend, defined using ", measure_txt, ",  over the time series,",
                                      " summarized over 1000 Monte-Carlo replications. Performance",
                                      " at sample size $n = 5000$ shown for ", all_est_txt, 
                                      " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                      " variable importance are displayed for each variable.", corr_txt,
                                      " \\label{tab:sim_bign_slope_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                     linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size) %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_slope_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
        
        # pick one estimator (SL)
        n_txt <- ifelse(this_dgm == 2, "$n \\in \\{250, 1000, 5000\\}$", "$n \\in \\{100, 250, 1000, 5000, 10000\\}$")
        summ_tab_sl <- summ_tib %>% 
          filter(Estimator == "SL") %>% 
          select(-outcome_type, -corr_between, -corr_within, -Estimator)
        knitr::kable(summ_tab_sl, 
                     format = "latex", booktabs = TRUE, escape = FALSE,
                     caption = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ",  and ",
                                      " the slope of the linear VIM trend over the time series,",
                                      " summarized over 1000 Monte-Carlo replications. Performance",
                                      " shown for each sample size ", n_txt, ",",
                                      " using the super learner.",
                                      " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                      " variable importance are displayed for each variable.", corr_txt,
                                      " \\label{tab:sim_sl_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                     linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size) %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_sl_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
        
        # separate average from slope
        summ_tib %>% 
          filter(grepl("average", `Summary`)) %>% 
          select(-outcome_type, -corr_between, -corr_within, -`Summary`) %>% 
          knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                       caption = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ",  over the time series,",
                                        " summarized over 1000 Monte-Carlo replications. Performance",
                                        " shown for each sample size ", n_txt, ",",
                                        " and ", all_est_txt,
                                        " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                        " variable importance are displayed for each variable.", corr_txt,
                                        " \\label{tab:sim_all_avg_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                       linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size, 
                                    latex_options = c("repeat_header"),
                                    repeat_header_text = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ",  over the time series,",
                                                                " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                    repeat_header_method = "replace") %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_avg_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
        summ_tib %>% 
          filter(grepl("slope", `Summary`)) %>% 
          select(-outcome_type, -corr_between, -corr_within, -`Summary`) %>% 
          knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                       caption = paste0("Performance of estimators of the linear VIM trend, defined using ", measure_txt, ",  over the time series,",
                                        " summarized over 1000 Monte-Carlo replications. Performance",
                                        " shown for each sample size ", n_txt, ",",
                                        " and ", all_est_txt, 
                                        " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                        " variable importance are displayed for each variable.", corr_txt,
                                        " \\label{tab:sim_all_slope_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                       linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size, 
                                    latex_options = c("repeat_header"),
                                    repeat_header_text = paste0("Performance of estimators of the linear VIM trend, defined using ", measure_txt, ",  over the time series,",
                                                                " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                    repeat_header_method = "replace") %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_slope_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
      }
    }
  }
  
}

# also print out separated by X to make the table less long
for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  for (j in 1:length(all_corrs)) {
    this_corr <- all_corrs[j]
    corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
    for (l in 1:length(all_measures)) {
      this_measure <- all_measures[l]
      this_nice_measure <- all_nice_measures[l]
      for (i in c(1, 2, 3, 8)) {
        this_summ_tib <- summ_tib %>% 
          filter(grepl(i, summ_tib$`VIM type: variable`), corr_within == this_corr,
                 Measure == this_nice_measure, dgm == this_dgm) %>% 
          select(-Measure, -dgm)
        all_est_txt <- paste0("all estimators: logistic regression (GLM),",
                              " lasso, random forests (RF), super learner (SL), and boosted trees (XGB)."
        )
        measure_txt <- ifelse(this_measure == "auc", "AUC", 
                              paste0(ifelse(this_measure == "ppv", "PPV", "sensitivity"), " at the 95th percentile of predicted risk"))
        if (nrow(this_summ_tib) == 0) {
          # do nothing
        } else {
          n_txt <- ifelse(this_dgm == 2, "$n \\in \\{250, 1000, 5000\\}$", "$n \\in \\{100, 250, 1000, 5000, 10000\\}$")
          this_summ_tib %>% 
            filter(grepl("average", `Summary`)) %>% 
            select(-outcome_type, -corr_between, -corr_within, -`Summary`) %>% 
            knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                         caption = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ", for variable ", i, " over the time series,",
                                          " summarized over 1000 Monte-Carlo replications. Performance",
                                          " shown for each sample size ", n_txt, ",",
                                          " and ", all_est_txt, 
                                          " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                          " variable importance are displayed.", corr_txt,
                                          " \\label{tab:sim_all_avg_", i, "_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                         linesep = "") %>% 
            kableExtra::kable_styling(font_size = font_size, 
                                      latex_options = c("repeat_header"),
                                      repeat_header_text = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ", for variable ", i, " over the time series,",
                                                                  " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                      repeat_header_method = "replace") %>% 
            kableExtra::save_kable(file = paste0(plots_dir, "summary_table_avg_", i, "_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
          this_summ_tib %>% 
            filter(grepl("slope", `Summary`)) %>% 
            select(-outcome_type, -corr_between, -corr_within, -`Summary`) %>% 
            knitr::kable(format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
                         caption = paste0("Performance of estimators of the linear VIM trend, defined using ", measure_txt, ", for variable ", i, " over the time series,",
                                          " summarized over 1000 Monte-Carlo replications. Performance",
                                          " shown for each sample size ", n_txt, ",",
                                          " and ", all_est_txt, 
                                          " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                          " variable importance are displayed.", corr_txt,
                                          " \\label{tab:sim_all_slope_", i, "_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                         linesep = "") %>% 
            kableExtra::kable_styling(font_size = font_size, 
                                      latex_options = c("repeat_header"),
                                      repeat_header_text = paste0("Performance of estimators of the linear VIM trend, defined using ", measure_txt, ", for variable ", i, " over the time series,",
                                                                  " summarized over 1000 Monte-Carlo replications. \\textit{(continued)}"),
                                      repeat_header_method = "replace") %>% 
            kableExtra::save_kable(file = paste0(plots_dir, "summary_table_slope_", i, "_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))  
        }
      }
    }
  }
  
}

summ_tab_bign <- summ_tib %>% 
  filter(`$n$` == 5000) %>% 
  select(-outcome_type, -corr_between, -`$n$`)


# summarize only for SL, GLM at n = 5000; each measure type
n_txt <- "$n = 5000$"
for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  variable_txt <- ifelse(this_dgm == 1, 
                         paste0("The true importance of variable 1",
                                " is nearly constant over time; variable 2 has increasing importance over time; variable 3",
                                " has decreasing importance over time; and variable 8 has zero importance at all time points."), 
                         paste0("The true importance of variable 1 is nearly constant over time;",
                                " variables 2 and 3 have increasing importance; and variable 8 has zero importance.")
  )
  for (i in 1:length(all_corrs)) {
    this_corr <- all_corrs[i]
    corr_txt <- ifelse(this_corr == 0, "", " Each feature has correlation 0.5 with the previous timepoint.")
    all_est_txt <- paste0("logistic regression (GLM) and super learner (SL).")
    for (j in 1:length(all_measures)) {
      this_measure <- all_measures[j]
      this_nice_measure <- all_nice_measures[j]
      measure_txt <- ifelse(this_measure == "auc", "AUC", 
                            paste0(ifelse(this_measure == "ppv", "PPV", "sensitivity"), " at the 95th percentile of predicted risk"))
      
      this_summ_tab_bign <- summ_tab_bign %>% 
        filter(corr_within == this_corr, Measure == this_nice_measure, dgm == this_dgm) %>% 
        select(-corr_within, -Measure, -dgm)
      
      if (nrow(this_summ_tab_bign) == 0) {
        # do nothing
      } else {
        big_n_summ <- this_summ_tab_bign %>% 
          filter(Estimator == "GLM" | Estimator == "SL") %>% 
          select(-Summary, -Estimator) %>% 
          knitr::kable(format = "latex", booktabs = TRUE,
                       caption = paste0("Performance of estimators of the average VIM, defined using ", measure_txt, ", and slope of the linear trend in VIM over the time series,",
                                        " summarized over 1000 Monte-Carlo replications. Performance",
                                        " at sample size ", n_txt, " shown for ", all_est_txt, 
                                        " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                        " variable importance are displayed for each variable.", variable_txt, corr_txt,
                                        " \\label{tab:sim_bign_all_", this_corr, "_", this_measure, "_dgm_", this_dgm, "}"),
                       linesep = "") %>% 
          kableExtra::kable_styling(font_size = font_size, 
                                    latex_options = c("scale_down")) %>% 
          kableExtra::pack_rows("Average (GLM)", 1, 8) %>% 
          kableExtra::pack_rows("Trend - slope (GLM)", 9, 16) %>% 
          kableExtra::pack_rows("Average (SL)", 17, 24) %>% 
          kableExtra::pack_rows("Trend - slope (SL)", 25, 32) 
        
        big_n_summ %>% 
          kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_all_", this_corr, "_", this_measure, "_dgm_", this_dgm, ".tex"))
      }
    }
  }
  
}

# summarize for all measure types, only SL, n = 5000
n_txt <- "$n = 5000$"
this_corr <- 0.5
all_measures_summ_bign <- summ_tab_bign %>% 
  filter(corr_within == 0.5, Estimator == "SL") %>% 
  select(-corr_within, -Estimator) %>% 
  arrange(Measure)

all_measure_txt <- paste0("PPV and sensitivity at the 95th percentile of predicted risk and AUC")

for (d in 1:length(all_dgms)) {
  this_dgm <- all_dgms[d]
  variable_txt <- ifelse(this_dgm == 1, 
                         paste0("The true importance of variable 1",
                                " is nearly constant over time; variable 2 has increasing importance over time; variable 3",
                                " has decreasing importance over time; and variable 8 has zero importance at all time points."), 
                         paste0("The true importance of variable 1 is nearly constant over time;",
                                " variables 2 and 3 have increasing importance; and variable 8 has zero importance.")
  )
  big_n_summ_all <- all_measures_summ_bign %>% 
    filter(dgm == this_dgm) %>% 
    select(-Measure, -Summary, -dgm) %>% 
    knitr::kable(format = "latex", booktabs = TRUE,
                 caption = paste0("\\small Performance of estimators of the average VIM, defined using ", all_measure_txt, ", and slope of the linear trend in VIM over the time series,",
                                  " summarized over 1000 Monte-Carlo replications. Performance",
                                  " at sample size ", n_txt, " shown for the super learner.",
                                  " Both add-in (compared to four covariates) and leave-out (compared to nine covariates)",
                                  " variable importance are displayed for each variable. ", variable_txt, corr_txt,
                                  " \\label{tab:sim_bign_all_measures_", this_corr, "_dgm_", this_dgm, "}"),
                 linesep = "") %>% 
    kableExtra::kable_styling(font_size = font_size - 2, 
                              latex_options = c("scale_down"))
  if (this_dgm == 1) {
      big_n_summ_all_final <- big_n_summ_all %>% 
        kableExtra::pack_rows("Average (AUC)", 1, 8) %>% 
        kableExtra::pack_rows("Trend - slope (AUC)", 9, 16)  
  } else {
    big_n_summ_all_final <- big_n_summ_all %>% 
      kableExtra::pack_rows("Average (AUC)", 1, 8) %>% 
      kableExtra::pack_rows("Trend - slope (AUC)", 9, 16) %>% 
      kableExtra::pack_rows("Average (PPV)", 17, 24) %>% 
      kableExtra::pack_rows("Trend - slope (PPV)", 25, 32) %>% 
      kableExtra::pack_rows("Average (Sensitivity)", 33, 40) %>% 
      kableExtra::pack_rows("Trend - slope (Sensitivity)", 41, 48)
  }
  
  
  big_n_summ_all_final %>% 
    kableExtra::save_kable(file = paste0(plots_dir, "summary_table_bign_all_measures_", this_corr, "_dgm_", this_dgm, ".tex"))
}

