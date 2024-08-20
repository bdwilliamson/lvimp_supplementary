# Running the numerical experiments for the `lvimp` paper

This file describes how to reproduce the experiments in ["Inference on summaries of a model-agnostic longitudinal variable importance trajectory with application to suicide prevention"](https://arxiv.org/abs/2311.01638) by Williamson, Moodie, Simon, Rossom, and Shortreed (_arXiv_, 2024+). All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. All analyses use the R package `lvimp` version 0.0.0.9000 and R package `vimp` version 2.3.3. 

The numerical experiments consist of two sections. First, we consider the properties of our proposal in a simple setting with no correlation among outcomes or covariates across timepoints. Second, we consider a more general setting with correlation among both the outcomes and covariates across time. However, both sets of results use identical code files, with varying tuning parameters. These files are:
* `utils.R`: generally-useful functions that are used across other files in the directory.
* `gen_data.R`: generate a dataset of a given sample size, number of timepoints, and correlation structure (among both outcomes and covariates). Other arguments include the number of predictors, the type of outcome (always set to `"binary"` in these experiments), and the true association between predictors and outcome and between confounding variables and the other variables.
* `investigate_cross_sectional_performance_once.R`: runs a single iteration of the simulation under a given set of parameters. Generates a dataset, runs a cross-validated Super Learner (with specified library of candidate learners) at each time point to estimate the outcome regression function, estimates cross-validated prediction performance and variable importance at each time point, and estimates longitudinal summaries of variable importance.
* `investigate_cross_sectional_performance.R`: runs the simulation a specified number of times for a given set of input parameters. Sets up varying parameters and fixed parameters.
* `compute_true_cross_sectional_values.R`: computes the true VIM values for the scenarios studied in the paper.
* `compile_cross_sectional_performance.R`: creates an analysis-ready dataset from the individual results files (from a given run of the simulation with fixed parameters).
* `plot_cross_sectional_performance.R`: creates plots and tables with results.

The file `create_figure_1.R` creates figure 1 from the paper, using simulated data.

## Results with no correlation across time

To reproduce results in the setting with no correlation, you can run the following code from a Windows command prompt:
```{cmd}
submit_cross_sectional_performance_uncorrelated.bat
```
This will run, sequentially, 1000 replications of the simulation in an uncorrelated setting for each sample size $n = 100, 250, 500, 1000, 5000, 10000$.

This Windows batch file runs `investigate_cross_sectional_performance.R` by setting the following arguments:
* `outcome-type`: binary (the outcome, here it is binary)
* `cor-between`: 0 (the correlation between predictors at a given time point)
* `cor-within`: 0 (the correlation within predictors and outcome across time)
* `n`: varies (the sample size)
* `p`: 10 (the number of predictors)
* `num-timepoints`: 4 (the number of timepoints)
* `nreps-total`: 1000 (the total desired number of simulation replications)
* `nreps-per-job`: 1000 (the number of simulation replications to run in each job)
* `parallel-cv`: 1 (runs cross-validation in parallel, the default)
* `run-leave-out`: 1 (runs leave-out VIMs)
* `run-add-in`: 1 (runs add-in VIMs)
* `simple-model`: 0 (run both simple and complex Super Learner algorithms)

## Results with correlation across time

To reproduce results in the setting with correlation, you can run the following code from a Windows command prompt:
```{cmd}
submit_cross_sectional_performance_correlated.bat
```
This will run, sequentially, 1000 replications of the simulation in an correlated setting for each sample size $n = 100, 250, 500, 1000, 5000, 10000$.

This Windows batch file runs `investigate_cross_sectional_performance.R` by setting the following arguments:
* `outcome-type`: binary (the outcome, here it is binary)
* `cor-between`: 0 (the correlation between predictors at a given time point)
* `cor-within`: 0.5 (the correlation within predictors and outcome across time)
* `n`: varies (the sample size)
* `p`: 10 (the number of predictors)
* `num-timepoints`: 4 (the number of timepoints)
* `nreps-total`: 1000 (the total desired number of simulation replications)
* `nreps-per-job`: 1000 (the number of simulation replications to run in each job)
* `parallel-cv`: 1 (runs cross-validation in parallel, the default)
* `run-leave-out`: 1 (runs leave-out VIMs)
* `run-add-in`: 1 (runs add-in VIMs)
* `simple-model`: 1 (run only simple algorithms in the Super Learner)

## Creating plots and tables

Once the simulations are finished, the plots and tables in the manuscript (and supplement) can be reproduced by running 

```{r}
compile_cross_sectional_performance.R
```

followed by

```{r}
plot_cross_sectional_performance.R
```