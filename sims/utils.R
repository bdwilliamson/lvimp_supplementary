# useful functions for longitudinal variable importance simulations

# measuring prediction performance for each learner in CV.SL -------------------
# get the CV-predictiveness for a single learner's predicted values
#' @param preds the fitted values
#' @param Y the outcome
#' @param full_y the observed outcome (from the entire dataset, for cross-fitted estimates)
#' @param scale what scale should the IPCW correction be applied on?
#'              Can help with numbers outside of (0, 1)
#'             ("identity" denotes the identity scale;
#'              "logit" means that AUC is transformed to the logit scale,
#'              correction is applied, and then back-transformed)
#' @param weights the inverse probability of censoring weights
#' @param C the indicator of being observed in phase 2 (1) or not (0)
#' @param Z a matrix of predictors observed on all participants in phase 1
#'          (can include the outcome)
#' @param type the type of predictiveness (e.g., auc, r-squared)
#' @param ... other arguments to measure_auc, but should include at least:
#'            a library of learners (using arg "SL.library")
#'            and may include control parameters for the super learner
#'            (e.g., cvControl = list(V = 5)
#'             for 5-fold cross-validated super learner)
one_predictiveness <- function(preds, Y, full_y = NULL, scale = "identity",
                    weights = rep(1, length(Y)), C = rep(1, length(Y)),
                    Z = NULL, type = "auc", ...) {
  if (type == "auc") {
    metric <- vimp::measure_auc
  } else {
    metric <- vimp::measure_r_squared
  }
  est_lst <- metric(
    fitted_values = preds, y = Y, full_y = full_y, C = C, Z = Z,
    ipc_weights = weights,
    ipc_fit_type = "SL", ...
  )
  list(est = est_lst$point_est, eif = est_lst$eif)
}

# get the cross-fitted CV-predictiveness for a single learner's predicted values
#' @param preds the fitted values
#' @param Y the outcome
#' @param folds the different cv folds that the learner was evaluated on
#' @param scale what scale should the IPCW correction be applied on?
#'              Can help with numbers outside of (0, 1)
#'             ("identity" denotes the identity scale;
#'              "logit" means that AUC is transformed to the logit scale,
#'              correction is applied, and then back-transformed)
#' @param weights the inverse probability of censoring weights
#' @param C the indicator of being observed in phase 2 (1) or not (0)
#' @param Z a matrix of predictors observed on all participants in phase 1
#'          (can include the outcome)
#' @param type the type of predictiveness (e.g., auc, r-squared)
#' @param ... other arguments to measure_auc, but should include at least:
#'            a library of learners (using arg "SL.library")
#'            and may include control parameters for the super learner
#'            (e.g., cvControl = list(V = 5)
#'             for 5-fold cross-validated super learner)
cv_predictiveness <- function(preds, Y, folds, scale = "identity",
                   weights = rep(1, length(Y)), C = rep(1, length(Y)),
                   Z = NULL, type = "auc", ...) {
  V <- length(folds)
  folds_numeric <- get_cv_sl_folds(folds)
  if (is.null(Z)) {
    folds_z <- folds_numeric
  } else {
    folds_z <- c(folds_numeric, sample(seq_len(V), nrow(Z) - length(folds_numeric),
                                       replace = TRUE
    ))
  }
  ests_eifs <- lapply(as.list(seq_len(V)), function(v) {
    one_predictiveness(
      preds = preds[folds_numeric == v], Y[folds_numeric == v],
      full_y = Y, scale = scale,
      weights = weights[folds_z == v], C = C[folds_z == v],
      Z = Z[folds_z == v, , drop = FALSE], type = type, ...
    )
  })
  est <- mean(unlist(lapply(ests_eifs, function(l) l$est)))
  var <- mean(unlist(lapply(ests_eifs, function(l) mean(l$eif ^ 2))))
  se <- sqrt(var / length(Y))
  ci <- vimp::vimp_ci(est, se, scale = scale, level = 0.95)
  return(list(est = est, se = se, ci = ci))
}

# get the folds from a CV.SL object, make them a vector
#' @param cv_sl_folds the CV.SL folds (a named list of row numbers)
#' @return a vector with the correct folds
get_cv_sl_folds <- function(cv_sl_folds) {
  folds_with_row_nums <- sapply(1:length(cv_sl_folds),
                                function(x) {
                                  list(
                                    row_nums = cv_sl_folds[[x]],
                                    fold = rep(x, length(cv_sl_folds[[x]]))
                                  )
                                },
                                simplify = FALSE
  )
  folds_df <- data.table::rbindlist(folds_with_row_nums)
  folds_df$fold[order(folds_df$row_nums)]
}
# given a vector, make a list of folds to pass to CV.SL
#' @param cv_folds a vector of CV folds
#' @return a named list of row numbers
make_cv_sl_folds <- function(cv_folds) {
  sorted_folds <- sort(unique(cv_folds))
  cv_sl_folds <- lapply(as.list(sorted_folds), function(v) {
    which(cv_folds == v)
  })
  names(cv_sl_folds) <- as.character(sorted_folds)
  return(cv_sl_folds)
}

# get the CV-predictiveness for an individual learner in the SL
#' @param sl_fit the fitted SL
#' @param col the column of interest (corresponds to a fitted algorithm)
#' @param scale scale for EIF estimation
#' @param weights IPC weights
#' @param C the censoring indicator
#' @param Z data observed on all participants
#' @param type the type of predictiveness (e.g., auc, r-squared)
#' @param ... other arguments to pass to Super Learner (for EIF estimation)
get_individual_predictiveness <- function(sl_fit, col, scale = "identity",
                               weights = rep(1, length(sl_fit$Y)),
                               C = rep(1, length(sl_fit$Y)), Z = NULL, type = "auc", ...) {
  if (any(is.na(sl_fit$library.predict[, col]))) {
    return(NULL)
  }
  alg_pred <- cv_predictiveness(
    preds = sl_fit$library.predict[, col], Y = sl_fit$Y,
    scale = scale,
    folds = sl_fit$folds, weights = weights,
    C = C, Z = Z, type = type, ...
  )
  # break apart algorithm and screen
  str <- colnames(sl_fit$library.predict)[col]
  alg <- gsub("(.*)_\\w+", "\\1", str)
  screen <- gsub("(.*)_{1}(.*)", "\\2", str)
  data.frame(
    Learner = alg, Screen = screen, est = alg_pred$est,
    se = alg_pred$se,
    ci_ll = alg_pred$ci[1], ci_ul = alg_pred$ci[2]
  )
}
# get the CV-predictiveness for all learners fit with SL
#' @param sl_fit the super learner fit object
#' @param scale what scale should the IPCW correction be applied on?
#'              Can help with numbers outside of (0, 1)
#'             ("identity" denotes the identity scale;
#'              "logit" means that AUC is transformed to the logit scale,
#'              correction is applied, and then back-transformed)
#' @param weights the inverse probability of censoring weights
#' @param C the indicator of being observed in phase 2 (1) or not (0)
#' @param Z a matrix of predictors observed on all participants in phase 1
#'          (can include the outcome)
#' @param type the type of predictiveness (e.g., auc, r-squared)
#' @param ... other arguments to measure_auc, but should include at least:
#'            a library of learners (using arg "SL.library")
#'            and may include control parameters for the super learner
#'            (e.g., cvControl = list(V = 5)
#'             for 5-fold cross-validated super learner)
get_all_predictiveness <- function(sl_fit, scale = "identity",
                         weights = rep(1, length(sl_fit$Y)),
                         C = rep(1, length(sl_fit$Y)),
                         Z = NULL, type = "auc", ...) {
  # get the CV-AUC of the SuperLearner predictions
  sl_pred <- cv_predictiveness(
    preds = sl_fit$SL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds,
    scale = scale, weights = weights, C = C, Z = Z, type = type, ...
  )
  out <- data.frame(
    Learner = "SL", Screen = "All", est = sl_pred$est,
    se = sl_pred$se, ci_ll = sl_pred$ci[1], ci_ul = sl_pred$ci[2]
  )
  
  # Get the CV-auc of the Discrete SuperLearner predictions
  discrete_sl_pred <- cv_predictiveness(
    preds = sl_fit$discreteSL.predict, Y = sl_fit$Y,
    folds = sl_fit$folds, scale = scale,
    weights = weights, C = C,
    Z = Z, type = type, ...
  )
  out <- rbind(out, data.frame(
    Learner = "Discrete SL", Screen = "All",
    est = discrete_sl_pred$est,
    se = discrete_sl_pred$se,
    ci_ll = discrete_sl_pred$ci[1],
    ci_ul = discrete_sl_pred$ci[2]
  ))
  
  # Get the cvauc of the individual learners in the library
  other_preds <- plyr::ldply(
    1:ncol(sl_fit$library.predict),
    function(x) {
      get_individual_predictiveness(
        sl_fit = sl_fit,
        col = x,
        scale = scale,
        weights = weights,
        C = C, Z = Z, type = type, ...
      )
    }
  )
  rbind(out, other_preds)
}
# extract best learner of a given type from CV.SL object -----------------------
#' @param learners a character vector specifying the different learners
get_unique_learners <- function(learners, unique_only = TRUE) {
  without_sl <- gsub("SL.", "", learners)
  without_tuning_parameters <- gsub("_.*", "", without_sl)
  if (unique_only) {
    return(unique(without_tuning_parameters))
  } else{
    return(without_tuning_parameters)
  }
}
#' Extract best learner of a given type from a CV.SL object
#' 
#' @param cv_sl the CV.SL object
#' @param procedure the type of procedure (e.g., "rf") for which to find the best learner
#' @param metric the prediction performance for each learner in the CV.SL object
#' 
#' @return the predictions based on the best learner from the given procedure type 
extract_best_learner <- function(cv_sl = NULL, procedure = "rf", prediction_performance = NULL) {
  all_learners <- get_unique_learners(prediction_performance$Learner, unique_only = FALSE)
  this_learner_perf <- prediction_performance[grepl(paste0("\\b", procedure, "\\b"), all_learners), ]
  # get the name of the learner with the best prediction performance
  best_learner <- this_learner_perf$Learner[which.max(this_learner_perf$est)]
  all_algs <-  gsub("(.*)_\\w+", "\\1", colnames(cv_sl$library.predict))
  # return predictions from that learner
  if (procedure != "SL") {
    return(list("learner" = best_learner, "preds" = cv_sl$library.predict[, best_learner == all_algs]))  
  } else {
    return(list("learner" = "SL", "preds" = cv_sl$SL.predict))
  }
  
}

# collapse list from simulation output -----------------------------------------
#' @param output_list a list of output for a given algorithm
#' @param algo the algorithm used for estimation
#' @param varsets the variable sets
#' @param vim_types the VIM types
#' @return a data.table with the correct output
collapse_sim_output <- function(output_list, algo = "glm", varsets = list(4:7), vim_types = "addi") {
  # extract point estimates, etc. of predictiveness over time for each variable set
  baseline_varset <- unlist(lapply(varsets, function(varset) length(setdiff(varset, 4:7)) == 0))
  all_varset <- unlist(lapply(varsets, function(varset) length(intersect(varset, 1:10)) == 10))
  vim_varsets <- varsets[!(baseline_varset | all_varset)]
  vim_vimtypes <- vim_types[!(baseline_varset | all_varset)]
  
  first_addi_indx <- which(unlist(vim_types) == "addi")[1]
  first_loco_indx <- which(unlist(vim_types) == "loco")[1]
  
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
  
  # cross_sectional_predictiveness <- rbind(
  #   data.table("varset" = "baseline", "designation" = paste0("predictiveness-", output_list[[1]]$timepoints),
  #              "est" = output_list[[1]]$predictiveness_reduced, 
  #              "se" = unlist(lapply(output_list[[1]]$vims, function(time_vim) time_vim$se_reduced)),
  #              "cil" = unlist(lapply(output_list[[1]]$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 1])),
  #              "ciu" = unlist(lapply(output_list[[1]]$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 2])),
  #              "p_value" = rep(NA, length(output_list[[1]]$timepoints))),
  #   data.table::rbindlist(lapply(output_list, function(l) {
  #     data.table("varset" = l$vims[[1]]$s, "designation" = paste0("predictiveness-", l$timepoints),
  #                "est" = l$predictiveness_full, "se" = unlist(lapply(l$vims, function(time_vim) time_vim$se_full)),
  #                "cil" = unlist(lapply(l$vims, function(time_vim) time_vim$predictiveness_ci_full[, 1])),
  #                "ciu" = unlist(lapply(l$vims, function(time_vim) time_vim$predictiveness_ci_full[, 2])),
  #                "p_value" = rep(NA, length(l$timepoints)))
  #   }))
  # )
  # extract point estimates, etc. of VIMs over time for each variable set
  cross_sectional_vims <- data.table::rbindlist(lapply(output_list, function(l) {
    data.table("varset" = l$vims[[1]]$s, "designation" = paste0("vim-", l$timepoints),
               "est" = l$vim, 
               "se" = unlist(lapply(l$vims, function(time_vim) time_vim$se)),
               "cil" = unlist(lapply(l$vims, function(time_vim) time_vim$ci[, 1])),
               "ciu" = unlist(lapply(l$vims, function(time_vim) time_vim$ci[, 2])),
               "p_value" = unlist(lapply(l$vims, function(time_vim) time_vim$p_value)))
  }))
  # extract summaries of predictiveness trajectories
  # summary_predictiveness <- rbind(
  #   data.table("varset" = "baseline", "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
  #              "est" = c(output_list[[1]]$average_reduced, output_list[[1]]$trend_reduced, output_list[[1]]$autc_reduced),
  #              "se" = c(output_list[[1]]$average_reduced_se, output_list[[1]]$trend_reduced_se, output_list[[1]]$autc_reduced_se),
  #              "cil" = c(output_list[[1]]$average_reduced_ci[, 1], output_list[[1]]$trend_reduced_ci[, 1], output_list[[1]]$autc_reduced_ci[, 1]),
  #              "ciu" = c(output_list[[1]]$average_reduced_ci[, 2], output_list[[1]]$trend_reduced_ci[, 2], output_list[[1]]$autc_reduced_ci[, 2]),
  #              "p_value" = c(output_list[[1]]$average_reduced_p_value, NA, output_list[[1]]$trend_reduced_p_value, output_list[[1]]$autc_reduced_p_value)),
  #   data.table::rbindlist(lapply(output_list, function(l) {
  #     data.table("varset" = l$vims[[1]]$s, "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
  #                "est" = c(l$average_full, l$trend_full, l$autc_full),
  #                "se" = c(l$average_full_se, l$trend_full_se, l$autc_full_se),
  #                "cil" = c(l$average_full_ci[, 1], l$trend_full_ci[, 1], l$autc_full_ci[, 1]),
  #                "ciu" = c(l$average_full_ci[, 2], l$trend_full_ci[, 2], l$autc_full_ci[, 2]),
  #                "p_value" = c(l$average_full_p_value, NA, l$trend_full_p_value, l$autc_full_p_value))
  #   }))
  # )
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
# extract the relevant output: point estimates, etc.
#' @param output the output for a given variable set
#' @param measure the measure to return: predictiveness, VIM, or summary
#' @param varset the variable set of interest
#' @param vim_type the type of variable importance (add-in or leave-out)
#' @return point estimate, SE, CI, and p-value
extract_predictiveness_output <- function(output, measure = "predictiveness", varset = c(4, 5, 6, 7),
                           vim_type = "loco") {
  is_predictiveness <- as.numeric(grepl("predictiveness", measure)) + 1
  desig_cross_section <- paste0(switch(is_predictiveness, "vim-", "predictiveness-"), output$timepoints)
  desig_trend <- paste0(switch(is_predictiveness, "vim-", "predictiveness-"), c("average", "trend-intercept", "trend-slope", "autc"))
  if (grepl("addi", vim_type)) { # need to keep track of predictiveness or VIM
    dt_cross_section <-   data.table("varset" = output$vims[[1]]$s, "designation" = desig_cross_section,
                                     "est" = output$predictiveness_full, "se" = unlist(lapply(output$vims, function(time_vim) time_vim$se_full)),
                                     "cil" = unlist(lapply(output$vims, function(time_vim) time_vim$predictiveness_ci_full[, 1])),
                                     "ciu" = unlist(lapply(output$vims, function(time_vim) time_vim$predictiveness_ci_full[, 2])),
                                     "p_value" = rep(NA, length(output$timepoints)))
    
    dt_summary <- data.table("varset" = output$vims[[1]]$s, "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
                             "est" = c(output$average_full, output$trend_full, output$autc_full),
                             "se" = c(output$average_full_se, output$trend_full_se, output$autc_full_se),
                             "cil" = c(output$average_full_ci[, 1], output$trend_full_ci[, 1], output$autc_full_ci[, 1]),
                             "ciu" = c(output$average_full_ci[, 2], output$trend_full_ci[, 2], output$autc_full_ci[, 2]),
                             "p_value" = c(output$average_full_p_value, NA, output$trend_full_p_value, output$autc_full_p_value))
  } else {
    dt_cross_section <- data.table("varset" = output$vims[[1]]$s, "designation" = desig_cross_section,
                                   "est" = output$predictiveness_reduced, "se" = unlist(lapply(output$vims, function(time_vim) time_vim$se_reduced)),
                                   "cil" = unlist(lapply(output$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 1])),
                                   "ciu" = unlist(lapply(output$vims, function(time_vim) time_vim$predictiveness_ci_reduced[, 2])),
                                   "p_value" = rep(NA, length(output$timepoints)))
    dt_summary <- data.table("varset" = output$vims[[1]]$s, "designation" = paste0("predictiveness-", c("average", "trend-intercept", "trend-slope", "autc")),
                             "est" = c(output$average_reduced, output$trend_reduced, output$autc_reduced),
                             "se" = c(output$average_reduced_se, output$trend_reduced_se, output$autc_reduced_se),
                             "cil" = c(output$average_reduced_ci[, 1], output$trend_reduced_ci[, 1], output$autc_reduced_ci[, 1]),
                             "ciu" = c(output$average_reduced_ci[, 2], output$trend_reduced_ci[, 2], output$autc_reduced_ci[, 2]),
                             "p_value" = c(output$average_reduced_p_value, NA, output$trend_reduced_p_value, output$autc_reduced_p_value))
  }
  return(list("cross_section" = dt_cross_section, "summary" = dt_summary))
}
# functions for nice plotting --------------------------------------------------
# turn variable set and designation into a nice plot suffix
#' @param varset the variable set, a character string
#' @param designation the designation, a character string
#' @return a string with the "nice" variable set and designation combination
nice_varset_designation <- function(varset = "baseline", designation = "predictiveness-1") {
  suffix <- ""
  if (grepl("baseline", varset)) {
    suffix <- paste0(suffix, "baseline variables")
  } else {
    suffix <- paste0(suffix, "variable ", varset, " (plus baseline)")
  }
  designation_vec <- unlist(strsplit(designation, "-", fixed = TRUE))
  suffix <- paste0(suffix, ", measure = ", designation_vec[1], ", timepoint = ", designation_vec[2])
  return(suffix)
}

# turn VIM prefix (e.g., "addi_" or "loco_") into a nice plot suffix
#' @param vim_descr the VIM of interest (addi_, for add-in; loco_, for leave-out) and the variable of interest
#' @param designation the designation (e.g., "predictiveness-1")
#' @return a string with the "nice" variable set and designation combination
nice_vim_designation <- function(vim_descr = "baseline", designation = "predictiveness-1") {
  suffix <- ""
  vim_descr_vec <- unlist(strsplit(vim_descr, "_", fixed = TRUE))
  designation_vec <- unlist(strsplit(designation, "-", fixed = TRUE))
  if (grepl("baseline", vim_descr)) {
    suffix <- paste0(suffix, "Predictiveness of baseline variables")
  } else {
    suffix <- paste0(suffix, ifelse(grepl("predictiveness", designation_vec[1]), "Predictiveness ", "VIM "), "of ")
    if (grepl("predictiveness", designation)) {
      suffix <- paste0(suffix, ifelse(grepl("addi", vim_descr), paste0("variable ", vim_descr_vec[2], " compared to baseline"),
                                      paste0("all variables besides variable ", vim_descr_vec[2])))
    } else {
      suffix <- paste0(suffix, "variable ", vim_descr_vec[2], ifelse(grepl("addi", vim_descr_vec[1]), " compared to baseline", " compared to all other variables"))  
    }
  }
  if (!is.na(suppressWarnings(as.numeric(tail(designation_vec, n = 1))))) {
    suffix <- paste0(suffix, ", timepoint = ", tail(designation_vec, n = 1))  
  } else {
    suffix <- paste0(suffix, ", ", ifelse(grepl("autc", designation), "AUTC of trajectory", ""),
                     ifelse(grepl("average", designation), "average over time", ""),
                     ifelse(grepl("trend-intercept", designation), "intercept of linear trend", ""),
                     ifelse(grepl("trend-slope", designation), "slope of linear trend", ""))
  }
  return(suffix)
}