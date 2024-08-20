# Analyzing the variable importance of predictors of suicide risk

This file describes the data analysis conducted in ["Inference on summaries of a model-agnostic longitudinal variable importance trajectory with application to suicide prevention"](https://arxiv.org/abs/2311.01638) by Williamson, Moodie, Simon, Rossom, and Shortreed (_arXiv_, 2024+). All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. All analyses use the R package `lvimp` version 0.0.0.9000 and R package `vimp` version 2.3.3. 

The datasets analyzed during this study are not publicly available because they contain detailed information from the electronic health records in the health systems participating in this study and are governed by Health Insurance Portability and Accountability Act (HIPAA). Data are, however, available from the authors upon reasonable request, with permission of all health systems involved and a fully executed data use agreement.

## Creating and describing the study cohort

As mentioned in the manuscript, our analysis dataset is a subset of a broader dataset. The file `01_create_cohort_study_analysis_dataset.R` creates this analysis dataset, by subsampling up to 6 visits per person spaced 90 days apart. The code in this file further selects variables of interest for the variable importance analysis and performs some quality checks.

The file `01_cohort_study_summary_statistics.R` creates Table S23, summarizing descriptive statistics of our sample.

## Estimating longitudinal VIMs of suicide risk

The file `02_estimate_prediction_performance.R` provides the code to replicate our analyses. It uses functions from `../sims/utils.R` and from `00_utils.R`. 

After reading in the dataset, we create the variable groups of interest for the VIM analysis. There are 19 variable groups (see Table S25 for details), which provide add-in or leave-out VIMs when the prediction performance of one group of variables is contrasted with another that has left out the variable(s) of interest.

We then estimate prediction performance for each variable set (using the Super Learner to estimate the prediction functions) at each of the six timepoints, and use these to estimate variable importance at each timepoint and summaries over time.

The file `02_estimate_prediction_performance_lasso_glm.R` is nearly identical, but uses a simplified Super Learner library consisting of only logistic regression and the lasso. From this analysis, we can inspect the individual algorithms in the Super Learner to see which variables are being selected by the lasso, for example, or the magnitude of regression coefficients in logistic regression models. This can help with diagnostics.

## Creating plots and tables

The file `03_plots_and_tables.R` contains code to create the plots and tables presented in the manuscript and supplement.