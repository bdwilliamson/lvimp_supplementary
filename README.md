# Supplementary materials for the `lvimp` paper

This repository contains code to reproduce the analyses in ["Inference on summaries of a model-agnostic longitudinal variable importance trajectory with application to suicide prevention"](https://arxiv.org/abs/2311.01638) by Williamson, Moodie, Simon, Rossom, and Shortreed (_arXiv_, 2024+). All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. All analyses use the R package `lvimp` version 0.0.0.9000 and R package `vimp` version 2.3.3. 

Further, all analyses were performed on multi-core virtual machines running the Windows operating system. There is code provided in each of the directories that will run the simulations from the Windows command prompt (`cmd`). If you use another system or batch job scheduler (e.g., a Linux cluster), these can serve as templates. You may also run the analyses locally, but they may take a large amount of time.

This README file provides an overview of the code available in the repository. We have separated our code into two directories based on the two main objectives of the manuscript. Each directory contains a separate README file that describes the code in that directory.

## The `sims` directory

This directory contains code to reproduce numerical experiments that describe the operating characteristics of our proposed method under varying data-generating mechanisms, in support of our main motivating data analysis. 

## The `data_analysis` directory

This directory contains code to reproduce the main data analysis, investigating the longitudinal variable importance of predictors of suicide risk, in particular the ninth item of the patient health questionnaire (PHQi9).

## Issues

If you encounter any bugs or have any specific questions about the analysis, please [file an issue](https://github.com/bdwilliamson/lvimp_supplementary/issues).