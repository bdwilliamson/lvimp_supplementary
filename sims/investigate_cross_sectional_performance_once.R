# Run the cross-sectional VIM simulation a single time
#' @param mc_id the monte-carlo id
#' @param n the sample size
#' @param p the number of features
#' @param outcome_type the type of outcome ("continuous" or "binary")
#' @param T the number of timepoints
#' @param beta_0 a list of matrices; the coefficients at each time point
#' @param confounder_beta the association between the confounders and the important variables
#' @param corr_between the between-variable correlation at a given time point
#' @param corr_within the within-variable correlation across time points
#' @param outcome_corr_type type of outcome correlation (defaults to "none", can be "ar" or "glmm")
#' @param learners a character vector of algorithms to pass to CV.SuperLearner
#' @param varsets the variable sets to estimate a conditional mean based upon
#' @param vim_types the types of variable importance ("addi" for add-in, "loco" for leave-out)
#' @param k_outer the number of folds for outer cross-fitting
#' @param k_inner the number of folds for inner cross-validation
#' @param parallel should we run CV.SuperLearner in parallel?
#' @param cl a cluster object (from parallel::makePSOCKcluster)
#' @return a tibble with results
investigate_cross_sectional_performance_once <- function(
  mc_id = 1, n = 100, p = 10, outcome_type = "binary", T = 4, 
  beta_0 = lapply(as.list(seq_len(args$num_timepoints)), function(t) {
    matrix(c(1, 1 + (t - 1) / 4, (-1) * (1 + exp((-1) * (t - 1)))^(-1) + 1, rep(0.25, 4), rep(0, args$p - 7)))
  }), confounder_beta = matrix(rep(0.05, 4)), corr_between = 0, corr_within = 0, 
  outcome_corr_type = "none",
  learners = "SL.glm", varsets = list(4:7), vim_types = "addi", k_outer = 10, k_inner = 10, 
  parallel = TRUE, cl = cl
) {
  # generate a dataset, in long format
  dataset <- gen_data(n = n, p = p, outcome_type = outcome_type, T = T,
                      beta_0 = beta_0, corr_between = corr_between, 
                      corr_within = corr_within, confounder_beta = confounder_beta)
  # estimate conditional means separately for each timepoint using CV.SuperLearner
  if (outcome_type == "binary") {
    sl_opts <- list("family" = "binomial", method = "method.CC_nloglik",
                    cvControl = list(V = k_outer, stratifyCV = TRUE),
                    innerCvControl = list(V = k_inner, stratifyCV = TRUE))
    measure_type <- "auc"
  } else {
    sl_opts <- list("family" = "gaussian", method = "method.CC_LS",
                    cvControl = list(V = k_outer),
                    innerCvControl = list(V = k_inner))
    measure_type <- "average_value"
  }
  unique_learners <- c(get_unique_learners(learners), "SL")
  cv_sls <- lapply(as.list(1:T), function(t) vector("list", length(varsets)))
  cv_folds <- vector("list", length = T)
  ss_folds <- lapply(as.list(1:T), function(i) vimp::make_folds(1:k_outer, V = 2))
  cv_pred_list <- lapply(as.list(1:T), function(t) vector("list", length(varsets)))
  cv_vim_list <- vector("list", length = T)
  for (i in 1:T) {
    these_timepoint_data <- dataset %>% filter(t == i)
    this_y <- these_timepoint_data$y
    this_x <- these_timepoint_data %>%
      select(-y, -t)
    cv_folds[[i]] <- vimp::make_folds(this_y, V = k_outer, stratified = TRUE)
    for (j in 1:length(varsets)) {
      if (!is.null(cv_folds[[i]])) {
        sl_opts$cvControl$validRows <- make_cv_sl_folds(cv_folds[[i]])
      }
      this_x_df <- this_x[, varsets[[j]]]
      if (parallel) {
        cv_sls[[i]][[j]] <- CV.SuperLearner(
          Y = this_y, X = as.matrix(this_x_df), SL.library = learners,
          family = sl_opts$family, method = sl_opts$method,
          cvControl = sl_opts$cvControl, innerCvControl = list(sl_opts$innerCvControl),
          parallel = cl
        ) 
      } else {
        cv_sls[[i]][[j]] <- CV.SuperLearner(
          Y = this_y, X = as.matrix(this_x_df), SL.library = learners,
          family = sl_opts$family, method = sl_opts$method,
          cvControl = sl_opts$cvControl, innerCvControl = list(sl_opts$innerCvControl),
          parallel = "seq"
        )
      }
    }
  }
  # estimate predictiveness, VIM for each group of variables at each timepoint
  # also estimate using each strategy; use best of each algorithm as assessed by CV
  confounder_varset_index <- which(unlist(lapply(varsets, function(x) {
    all_equal <- all.equal(x, 4:7)
    if (is.character(all_equal)) {
      return(FALSE)
    } else if (!all_equal) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  })))
  all_varset_index <- which(unlist(lapply(varsets, function(x) {
    all_equal <- all.equal(x, 1:10)
    if (is.character(all_equal)) {
      return(FALSE)
    } else if (!all_equal) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  })))
  non_confounder_all_varsets <- seq_len(length(varsets))[-c(confounder_varset_index, all_varset_index)]
  for (i in 1:T) {
    these_timepoint_data <- dataset %>% filter(t == i)
    this_y <- these_timepoint_data$y
    this_x <- these_timepoint_data %>%
      select(-y, -t)
    cv_folds_vec <- cv_folds[[i]]
    cv_vim_list[[i]] <- vector("list", length = length(non_confounder_all_varsets))
    # get confounder-only, all preds for each algorithm type
    # estimate CV prediction performance
    if (length(confounder_varset_index) == 0) {
      cv_pred_perf_confounders <- NA
      confounder_preds_list <- vector("list", length = length(unique_learners))
    } else {
      cv_pred_perf_confounders <- get_all_predictiveness(sl_fit = cv_sls[[i]][[confounder_varset_index]],
                                                         type = measure_type)
      confounder_preds_list <- vector("list", length = length(unique_learners))
      # get best algorithm within each class
      for (k in 1:length(unique_learners)) {
        this_best_learner <- extract_best_learner(
          cv_sl = cv_sls[[i]][[confounder_varset_index]], procedure = unique_learners[k],
          prediction_performance = cv_pred_perf_confounders
        )
        confounder_preds_list[[k]] <- this_best_learner$preds
        names(confounder_preds_list)[k] <- this_best_learner$learner
      }
    }
    if (length(all_varset_index) == 0) {
      cv_pred_perf_all <- NA
      all_preds_list <- vector("list", length = length(unique_learners))
    } else {
      cv_pred_perf_all <- get_all_predictiveness(sl_fit = cv_sls[[i]][[all_varset_index]],
                                                 type = measure_type)
      all_preds_list <- vector("list", length = length(unique_learners))
      # get best algorithm within each class
      for (k in 1:length(unique_learners)) {
        this_best_learner <- extract_best_learner(
          cv_sl = cv_sls[[i]][[all_varset_index]], procedure = unique_learners[k],
          prediction_performance = cv_pred_perf_all
        )
        all_preds_list[[k]] <- this_best_learner$preds
        names(all_preds_list)[k] <- this_best_learner$learner
      }
    }
    for (j in non_confounder_all_varsets) {
      # get the predictions from the best algorithm within each class
      this_preds_list <- vector("list", length = length(unique_learners))
      this_cv_pred_perf <- get_all_predictiveness(sl_fit = cv_sls[[i]][[j]],
                                                  type = measure_type)
      for (k in 1:length(unique_learners)) {
        this_best_learner <- extract_best_learner(
          cv_sl = cv_sls[[i]][[j]], procedure = unique_learners[k],
          prediction_performance = this_cv_pred_perf
        )
        this_preds_list[[k]] <- this_best_learner$preds
        names(this_preds_list)[k] <- this_best_learner$learner
      }
      # get estimated variable importance using the best algorithm in each class
      if (grepl("addi", vim_types[[j]])) {
        full_preds <- this_preds_list
        reduced_preds <- confounder_preds_list
      } else {
        full_preds <- all_preds_list
        reduced_preds <- this_preds_list
      }
      cv_vim_list[[i]][[which(j == non_confounder_all_varsets)]] <- lapply(as.list(1:length(unique_learners)), function(k) {
        this_vim <- vimp::cv_vim(Y = this_y, X = this_x,
                                 cross_fitted_f1 = full_preds[[k]],
                                 cross_fitted_f2 = reduced_preds[[k]],
                                 indx = varsets[[j]][!(varsets[[j]] %in% c(4, 5, 6, 7))],
                                 type = measure_type,
                                 cross_fitting_folds = cv_folds_vec,
                                 sample_splitting_folds = ss_folds[[i]],
                                 run_regression = FALSE,
                                 alpha = 0.05, V = k_outer / 2, 
                                 na.rm = TRUE)
        this_alg_full <- names(this_preds_list)[k]
        this_alg_reduced <- names(confounder_preds_list)[k]
        tibble("full" = this_alg_full, "reduced" = this_alg_reduced, "vim" = this_vim)
      })
    }
  }
  # reorganize the list to be in the order algorithm, variable set, timepoint
  reordered_cv_vim_list <- vector("list", length = length(unique_learners))
  for (k in 1:length(unique_learners)) {
    reordered_cv_vim_list[[k]] <- vector("list", length = length(non_confounder_all_varsets))
    for (j in 1:(length(cv_vim_list[[1]]))) {
      reordered_cv_vim_list[[k]][[j]] <- lapply(cv_vim_list, function(l) l[[j]][[k]])
    }
  }
  # estimate summaries of predictiveness, VIM for each algorithm and group of variables
  lvim_list <- vector("list", length = length(unique_learners))
  for (k in 1:length(unique_learners)) {
    lvim_list[[k]] <- vector("list", length = length(non_confounder_all_varsets))
    for (j in 1:length(reordered_cv_vim_list[[k]])) {
      this_vim_list <- lapply(reordered_cv_vim_list[[k]][[j]], function(l) l$vim)
      this_full_alg <- unlist(lapply(reordered_cv_vim_list[[k]][[j]], function(l) l$full[1]))
      this_reduced_alg <- unlist(lapply(reordered_cv_vim_list[[k]][[j]], function(l) l$reduced[1]))
      lvim_obj <- lvimp::lvim(this_vim_list, timepoints = 1:4)
      lvim_list[[k]][[j]] <- lvimp::lvim_average(lvim_obj, indices = 1:4)
      lvim_list[[k]][[j]] <- lvimp::lvim_trend(lvim_list[[k]][[j]], indices = 1:4)
      lvim_list[[k]][[j]] <- lvimp::lvim_autc(lvim_list[[k]][[j]], indices = 1:4)
    }
  }
  # return output!
  output <- cbind("mc_id" = mc_id, "n" = n, "p" = p, "outcome_type" = outcome_type, 
                  "corr_between" = corr_between, "corr_within" = corr_within,
                  data.table::rbindlist(
                    lapply(as.list(1:length(unique_learners)), function(i) {
                      collapse_sim_output(output_list = lvim_list[[i]], algo = unique_learners[i], varsets = varsets,
                                          vim_types = vim_types)
                    })
                  ))
  return(output)
}