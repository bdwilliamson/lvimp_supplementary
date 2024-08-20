# useful functions for the data analysis

# for getting predictions and prediction performance ---------------------------
# get predictions from the best learner of each class and the discrete and ensemble super learners
# @param cv_sl the CV.SuperLearner object
# @param learners the learners of interest
# @param perf_obj the prediction performance object
# @return the predictions corresponding to the best-in-class learner
get_all_learner_preds <- function(cv_sl = NULL, learners = "SL", perf_obj = NULL) {
  all_preds_list <- lapply(as.list(1:length(learners)), function(k) {
    best_learner <- extract_best_learner(
      cv_sl = cv_sl, procedure = learners[k], prediction_performance = perf_obj
    )
    best_learner$preds
  })
  names(all_preds_list) <- learners
  return(all_preds_list)
}

# get VIM based on the best learner in each class, discrete SL, and/or SL
# @param y the outcome
# @param x the covariates
# @param full_preds the predictions based on the larger set of covariates
# @param reduced_preds the predictions based on the smaller set of covariates
# @param var_set the variable set (index set)
# @param measure_type the type of VIM (e.g., AUC)
# @param cv_folds the cross-fitting folds
# @param ss_folds the sample-splitting folds
# @param alpha the type I error rate
# @param K the number of CV folds (for VIM estimation)
# @param complete_obs for longitudinal VIMs: the observations that are available at *all* timepoints
# @return a list with estimated VIMs for each algorithm
get_all_vims <- function(y = NULL, x = NULL, full_preds = NULL,
                         reduced_preds = NULL, var_set = 1, 
                         measure_type = "auc", cv_folds = rep(1, length(y)),
                         ss_folds = rep(1, length(y)), alpha = 0.05, K = 1,
                         complete_obs = rep(1, length(y))) {
  cv_vim_list <- lapply(as.list(1:length(full_preds)), function(k) {
    if (all(is.na(full_preds[[k]])) | all(is.na(reduced_preds[[k]]))) {
      this_vim <- list("s" = var_set, "est" = NA, "naive" = NA,
                       "eif" = NA, "eif_full" = NA, "eif_redu" = NA,
                       "ci" = matrix(c(NA, NA), nrow = 1), "se" = NA, "predictiveness_full" = NA,
                       "predictiveness_reduced" = NA, "predictiveness_ci_full" = matrix(c(NA, NA), nrow = 1),
                       "predictiveness_ci_reduced" = matrix(c(NA, NA), nrow = 1), "se_full" = NA,
                       "se_reduced" = NA, "test" = NA, "p_value" = NA,
                       "mat" = tibble::tibble(
                         "s" = var_set, "est" = NA, "se" = NA, "cil" = NA, "ciu" = NA,
                         "test" = NA, "p_value" = NA
                       ), "scale" = "identity", "alpha" = 0.05)
      class(this_vim) <- c("vim", "list")
    } else {
      this_vim <- vimp::cv_vim(Y = y, X = x, cross_fitted_f1 = full_preds[[k]],
                               cross_fitted_f2 = reduced_preds[[k]],
                               indx = var_set, type = measure_type,
                               cross_fitting_folds = cv_folds, ss_folds = ss_folds,
                               run_regression = FALSE, alpha = alpha,
                               V = K, na.rm = TRUE)
      fixed_se <- vimp::vimp_se(
        eif_full = lapply(this_vim$all_eifs_full, function(eif) eif[complete_obs]),
        eif_reduced = lapply(this_vim$all_eifs_redu, function(eif) eif[complete_obs]),
        cross_fit = TRUE, sample_split = TRUE, na.rm = FALSE
      )
      this_vim$se <- fixed_se
      var_full <- mean(
        unlist(lapply(this_vim$all_eifs_full, function(eif) mean(eif[complete_obs] ^ 2)))
      )
      this_vim$se_full <- sqrt(var_full / sum(complete_obs))
      var_reduced <- mean(
        unlist(lapply(this_vim$all_eifs_redu, function(eif) mean(eif[complete_obs] ^ 2)))
      )
      this_vim$se_reduced <- sqrt(var_reduced / sum(complete_obs))
      this_vim$ci <- vimp::vimp_ci(est = this_vim$est, se = this_vim$se, truncate = FALSE)
      this_vim$predictiveness_ci_full <- vimp::vimp_ci(est = this_vim$predictiveness_full,
                                                       se = this_vim$se_full, truncate = FALSE)
      this_vim$predictiveness_ci_reduced <- vimp::vimp_ci(est = this_vim$predictiveness_reduced,
                                                          se = this_vim$se_reduced, truncate = FALSE)
    }
    this_alg_full <- names(full_preds)[k]
    this_alg_reduced <- names(reduced_preds)[k]
    tibble("full" = this_alg_full, "reduced" = this_alg_reduced, "vim" = this_vim)
  })
  return(cv_vim_list)
}

# get all longitudinal vim summaries
# @param cv_vims a list of CV VIM estimates
# @param num_timepoints the number of timepoints
# @return longitudinal VIM summaries!
get_all_lvims <- function(cv_vims = NULL, num_timepoints = 4) {
  K <- length(cv_vims)
  num_varsets <- length(cv_vims[[1]])
  lvim_list <- vector("list", length = K)
  # outer loop over algorithms
  for (k in seq_len(K)) {
    lvim_list[[k]] <- vector("list", length = num_varsets)
    # inner loop over variable sets
    for (j in seq_len(num_varsets)) {
      this_vim_list <- lapply(cv_vims[[k]][[j]], function(l) l$vim)
      this_full_alg <- unlist(lapply(cv_vims[[k]][[j]], function(l) l$full[1]))
      this_reduced_alg <- unlist(lapply(cv_vims[[k]][[j]], function(l) l$reduced[1]))
      lvim_obj <- lvimp::lvim(this_vim_list, timepoints = 1:num_timepoints)
      lvim_list[[k]][[j]] <- lvimp::lvim_average(lvim_obj, indices = 1:num_timepoints)
      lvim_list[[k]][[j]] <- lvimp::lvim_trend(lvim_list[[k]][[j]], indices = 1:num_timepoints)
      lvim_list[[k]][[j]] <- lvimp::lvim_autc(lvim_list[[k]][[j]], indices = 1:num_timepoints)    
    }
  }
  return(lvim_list)
}

# collapse list output
#' @param output_list a list of output for a given algorithm
#' @param algo the algorithm used for estimation
#' @param varsets the variable sets
#' @param vim_types the VIM types
#' @return a data.table with the correct output
collapse_output <- function(output_list, algo = "glm", varsets = list(4:7), vim_types = "addi",
                                baseline_vars = 4:7, all_vars = 1:10) {
  # extract point estimates, etc. of predictiveness over time for each variable set
  baseline_varset <- unlist(lapply(varsets, function(varset) length(setdiff(varset, baseline_vars)) == 0))
  # all_varset <- unlist(lapply(varsets, function(varset) length(intersect(varset, all_vars)) == length(all_vars)))
  all_varset <- rep(FALSE, length(varsets))
  # vim_varsets <- varsets[!(baseline_varset | all_varset)]
  vim_varsets <- varsets[!baseline_varset]
  # vim_vimtypes <- vim_types[!(baseline_varset | all_varset)]
  vim_vimtypes <- vim_types[!baseline_varset]
  
  first_addi_indx <- which(unlist(vim_types) == "addi")[1]

  predictiveness_list <- lapply(as.list(1:length(output_list)), function(v) {
    extract_predictiveness_output(output = output_list[[v]], measure = "predictiveness",
                                  varset = vim_varsets[[v]], vim_type = vim_vimtypes[[v]])
  })
  cross_sectional_predictiveness <- data.table::rbindlist(lapply(predictiveness_list, function(x) x$cross_section))
  summary_predictiveness <- data.table::rbindlist(lapply(predictiveness_list, function(x) x$summary))
  if (any(baseline_varset)) {
    cross_sectional_predictiveness <- rbind(
      data.table("varset" = "baseline", "designation" = paste0("predictiveness-", output_list[[1]]$timepoints),
                 "est" = output_list[[first_addi_indx]]$predictiveness_reduced, 
                 "se" = unlist(lapply(output_list[[first_addi_indx]]$vims, function(time_vim) time_vim$se_reduced)),
                 "cil" = unlist(lapply(output_list[[first_addi_indx]]$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 1])),
                 "ciu" = unlist(lapply(output_list[[first_addi_indx]]$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 2])),
                 "p_value" = rep(NA, length(output_list[[first_addi_indx]]$timepoints))),
      cross_sectional_predictiveness
    )
    summary_predictiveness <- rbind(
      data.table("varset" = "baseline", "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
                 "est" = c(output_list[[first_addi_indx]]$average_reduced, output_list[[first_addi_indx]]$trend_reduced, output_list[[first_addi_indx]]$autc_reduced),
                 "se" = c(output_list[[first_addi_indx]]$average_reduced_se, output_list[[first_addi_indx]]$trend_reduced_se, output_list[[first_addi_indx]]$autc_reduced_se),
                 "cil" = c(output_list[[first_addi_indx]]$average_reduced_ci[, 1], output_list[[first_addi_indx]]$trend_reduced_ci[, 1], output_list[[first_addi_indx]]$autc_reduced_ci[, 1]),
                 "ciu" = c(output_list[[first_addi_indx]]$average_reduced_ci[, 2], output_list[[first_addi_indx]]$trend_reduced_ci[, 2], output_list[[first_addi_indx]]$autc_reduced_ci[, 2]),
                 "p_value" = c(output_list[[first_addi_indx]]$average_reduced_p_value, NA, output_list[[first_addi_indx]]$trend_reduced_p_value, output_list[[first_addi_indx]]$autc_reduced_p_value)),
      summary_predictiveness
    )
  }
  if (any(all_varset)) {
    cross_sectional_predictiveness <- rbind(
      cross_sectional_predictiveness,
      data.table("varset" = "all", "designation" = paste0("predictiveness-", output_list[[1]]$timepoints),
                 "est" = output_list[[first_loco_indx]]$predictiveness_full, 
                 "se" = unlist(lapply(output_list[[first_loco_indx]]$vims, function(time_vim) time_vim$se_full)),
                 "cil" = unlist(lapply(output_list[[first_loco_indx]]$vims, function(time_vim) time_vim$predictiveness_ci_full[, 1])),
                 "ciu" = unlist(lapply(output_list[[first_loco_indx]]$vims, function(time_vim) time_vim$predictiveness_ci_full[, 2])),
                 "p_value" = rep(NA, length(output_list[[first_loco_indx]]$timepoints)))
    )
    summary_predictiveness <- rbind(
      summary_predictiveness,
      data.table("varset" = "all", "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
                 "est" = c(output_list[[first_loco_indx]]$average_reduced, output_list[[first_loco_indx]]$trend_reduced, output_list[[first_loco_indx]]$autc_reduced),
                 "se" = c(output_list[[first_loco_indx]]$average_reduced_se, output_list[[first_loco_indx]]$trend_reduced_se, output_list[[first_loco_indx]]$autc_reduced_se),
                 "cil" = c(output_list[[first_loco_indx]]$average_reduced_ci[, 1], output_list[[first_loco_indx]]$trend_reduced_ci[, 1], output_list[[first_loco_indx]]$autc_reduced_ci[, 1]),
                 "ciu" = c(output_list[[first_loco_indx]]$average_reduced_ci[, 2], output_list[[first_loco_indx]]$trend_reduced_ci[, 2], output_list[[first_loco_indx]]$autc_reduced_ci[, 2]),
                 "p_value" = c(output_list[[first_loco_indx]]$average_reduced_p_value, NA, output_list[[first_loco_indx]]$trend_reduced_p_value, output_list[[first_loco_indx]]$autc_reduced_p_value))
    )
  }
  # extract point estimates, etc. of VIMs over time for each variable set
  cross_sectional_vims <- data.table::rbindlist(lapply(output_list, function(l) {
    data.table("varset" = l$vims[[1]]$s, "designation" = paste0("vim-", l$timepoints),
               "est" = l$vim, 
               "se" = unlist(lapply(l$vims, function(time_vim) time_vim$se)),
               "cil" = unlist(lapply(l$vims, function(time_vim) time_vim$ci[, 1])),
               "ciu" = unlist(lapply(l$vims, function(time_vim) time_vim$ci[, 2])),
               "p_value" = unlist(lapply(l$vims, function(time_vim) time_vim$p_value)))
  }))
  # extract summaries of VIM trajectories
  summary_vims <- data.table::rbindlist(lapply(output_list, function(l) {
    data.table("varset" = l$vims[[1]]$s, "designation" = paste0("vim-", c("average", "trend-intercept", "trend-slope", "autc")),
               "est" = c(l$average_vim, l$trend_vim, l$autc_vim),
               "se" = c(l$average_vim_se, l$trend_vim_se, l$autc_vim_se),
               "cil" = c(l$average_vim_ci[, 1], l$trend_vim_ci[, 1], l$autc_vim_ci[, 1]),
               "ciu" = c(l$average_vim_ci[, 2], l$trend_vim_ci[, 2], l$autc_vim_ci[, 2]),
               "p_value" = c(l$average_vim_p_value, NA, l$trend_vim_p_value, l$autc_vim_p_value))
  }))
  # return
  output <- cbind("algo" = algo, rbind(cross_sectional_predictiveness, cross_sectional_vims, 
                                       summary_predictiveness, summary_vims))
  return(output)
}

# for making plots and tables --------------------------------------------------
# get a nice and numeric variable set from a list of indices
# @param varset the variable index set
# @return the number corresponding to varset
get_varset_num <- function(varset = NULL) {
  case_when(
    is.null(varset) | varset == "baseline" ~ 1,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216" ~ 2,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130" ~ 3,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,66,67,68,69,70,71,72,73,74,75,76,77,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,121,122,123,124,125,126,127,128,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179" ~ 4,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,5" ~ 5,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,180" ~ 6,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,66,67,68,69,70,71,72,73,74,75,76,77,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,121,122,123,124,125,126,127,128,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,5" ~ 7,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,66,67,68,69,70,71,72,73,74,75,76,77,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,121,122,123,124,125,126,127,128,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180" ~ 8,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,5,180" ~ 9,
    varset == "1,2,3,4,6,7,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,66,67,68,69,70,71,72,73,74,75,76,77,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,121,122,123,124,125,126,127,128,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,5,180" ~ 10,
    varset == "180" ~ 11,
    varset == "3" ~ 12,
    varset == "3,184,185,186" ~ 13,
    varset == "3,184,185,186,187,188,189,190,191,192,193" ~ 14,
    varset == "3,180" ~ 15,
    varset == "3,184,185,186,180" ~ 16,
    varset == "3,184,185,186,187,188,189,190,191,192,193,180" ~ 17,
    varset == "3,184,185,186,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130" ~ 18,
    varset == "3,184,185,186,33,34,35,36,37,38,63,64,65,78,79,80,99,100,101,102,119,120,129,130,180" ~ 19
  )
}
# @param varset_num the variable set number
# @param varsets the vector of nice varsets
# @return the nice variable set
get_nice_varset <- function(varset_num, varsets) {
  case_when(
    varset_num == 1 ~ varsets[1],
    varset_num == 2 ~ varsets[2],
    varset_num == 3 ~ varsets[3],
    varset_num == 4 ~ varsets[4],
    varset_num == 5 ~ varsets[5],
    varset_num == 6 ~ varsets[6],
    varset_num == 7 ~ varsets[7],
    varset_num == 8 ~ varsets[8],
    varset_num == 9 ~ varsets[9],
    varset_num == 10 ~ varsets[10],
    varset_num == 11 ~ varsets[11],
    varset_num == 12 ~ varsets[12],
    varset_num == 13 ~ varsets[13],
    varset_num == 14 ~ varsets[14],
    varset_num == 15 ~ varsets[15],
    varset_num == 16 ~ varsets[16],
    varset_num == 17 ~ varsets[17],
    varset_num == 18 ~ varsets[18],
    varset_num == 19 ~ varsets[19]
  )
}

# @param output the output dataset
# @param varset the variable set(s) of interest
# @param est_type the designation (predictiveness or vim)
# @param use_paper_numbering should we use paper numbering or original numbering?
# @return a ggplot object with the trajectory plotted over time
plot_trajectory <- function(output, varset = 1, est_type = "predictiveness",
                            use_paper_numbering = TRUE) {
  if (use_paper_numbering) {
    this_output <- output %>% 
      filter(paper_varset_num %in% varset, grepl(est_type, designation) & !grepl("-[average|trend|autc]", designation)) %>% 
      mutate(timepoint = gsub(paste0(est_type, "-"), "", designation))   
  } else {
    this_output <- output %>% 
      filter(varset_num %in% varset, grepl(est_type, designation) & !grepl("-[average|trend|autc]", designation)) %>% 
      mutate(timepoint = gsub(paste0(est_type, "-"), "", designation)) 
  }
  if (grepl("vim", est_type)) {
    this_output <- this_output %>% 
      mutate(est = pmax(0, est), cil = pmax(0, cil))
  }
  if (length(varset) > 1) {
    plot_init <- this_output %>% 
      ggplot(aes(x = timepoint, y = est, shape = algo_fct, color = varset_fct))
  } else {
    plot_init <- this_output %>% 
      ggplot(aes(x = timepoint, y = est, shape = algo_fct))
  } 
  y_lim <- switch(as.numeric(grepl("vim", est_type)) + 1, c(0.35, 1), c(0, 0.5))
  this_plot <- plot_init +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.5)) +
    scale_color_viridis_d(begin = 0, end = 0.75) +
    ggtitle(paste0(toupper(est_type), " TRAJECTORY")) +
    labs(shape = "Algorithm", x = "Time point", y = "AUC") +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(shape = guide_legend(nrow = 2)) +
    ylim(y_lim)
  if (length(varset) > 1) {
    this_plot <- this_plot + 
      labs(color = "Variable set") +
      guides(color = guide_legend(nrow = 2))
  }
  return(this_plot)
}

# @param output the output dataset
# @param varset the variable set of interest
# @param est_type the designation, vim or predictiveness
# @param trajectory whether to look at the trajectory (true) or not (false)
# @param min_zero is the minimum value zero? 
# @return a table with the longitudinal summaries
create_overall_table <- function(output, varset = 1, est_type = "predictiveness",
                                  trajectory = TRUE, digits = 3,
                                 min_zero = TRUE) {
  if (trajectory) {
    this_output <- output %>% 
      filter(varset_num %in% varset, grepl(est_type, designation) & !grepl("-[average|trend|autc]", designation)) %>% 
      mutate(measure = gsub(paste0(est_type, "-"), "", designation))
  } else {
    this_output <- output %>% 
      filter(varset_num %in% varset, grepl(est_type, designation) & grepl("-[average|trend|autc]", designation)) %>% 
      mutate(measure = gsub(paste0(est_type, "-"), "", designation))
  }
  # round point estimate, SE, CI (and make one column), return only necessary columns
  if (min_zero) {
    nice_output <- this_output %>% 
      mutate(est = ifelse(est < 0 & trajectory, 0, round(est, digits)),
             cil = ifelse(cil < 0 & trajectory, 0, round(cil, digits)))
  } else {
    nice_output <- this_output %>% 
      mutate(est = round(est, digits), cil = round(cil, digits))
  }
  nice_output <- nice_output %>% 
    mutate(ciu = round(ciu, digits),
           ci = paste0("[", cil, ", ", ciu, "]"),
           p_value = ifelse(p_value < 0.001, "< 0.001", as.character(round(p_value, digits)))) %>% 
    arrange(measure, varset_fct, algo_fct) %>% 
    select(measure, nice_varset, algo, est, se, ci, p_value) %>% 
    rename(Algorithm = algo, `Point Estimate` = est, `SE` = se, `95% CI` = ci,
           `p-value` = p_value, `Variable set` = nice_varset)
  if (grepl("predictiveness", est_type)) {
    ret <- nice_output %>% 
      select(-`p-value`)
  } else {
    ret <- nice_output
  }
  if (trajectory) {
    ret <- ret %>% rename(Timepoint = measure)
  } else {
    ret <- ret %>% 
      mutate(measure = case_when(
        measure == "autc" ~ "AUTC",
        measure == "trend-intercept" ~ "Intercept (linear trend)",
        measure == "trend-slope" ~ "Slope (linear trend)",
        measure == "average" ~ "Average"
      )) %>% 
      rename(`Summary measure` = measure)
  }
  return(ret)
}

# @param output the output dataset
# @param varset the variable set of interest
# @param comparison 
# @return a ggplot object with the trajectory plotted over time
create_comparison_table <- function(output, varset = 1, est_type = "predictiveness", 
                                    comparison = "lag") {
  
}

