# compile output from all simulations

# load required packages -------------------------------------------------------
library("data.table")
library("here")
library("dplyr")

# read in results --------------------------------------------------------------
results_dir <- here::here("..", "results", "sims")
results_dir <- "H:/Papers/longitudinal_vim/results/sims/"
# note results are based on compiled results (i.e., each combination has 1000 replicates)
all_files <- list.files(results_dir, pattern = "output_[binary|continuous]*")
all_files <- all_files[!grepl("interim", all_files)]
all_output_lst <- lapply(as.list(all_files), function(file) {
  readRDS(paste0(results_dir, file))
})
all_output <- data.table::rbindlist(all_output_lst, fill = TRUE)
# fix varset
all_output_2 <- all_output %>% 
  mutate(numeric_varset = suppressWarnings(as.numeric(varset))) %>% 
  mutate(varset = case_when(
    varset == "baseline" ~ "baseline",
    varset == "all" ~ "all",
    !is.na(numeric_varset) ~ paste0("addi_", numeric_varset),
    varset == "2,3,8,9,10" ~ "loco_1",
    varset == "1,3,8,9,10" ~ "loco_2",
    varset == "1,2,8,9,10" ~ "loco_3",
    varset == "1,2,3,9,10" ~ "loco_8",
    varset == "1,2,3,8,10" ~ "loco_9",
    varset == "1,2,3,8,9" ~ "loco_10"
  )) %>% 
  select(-numeric_varset)
# attach on the true values
all_truth_files <- list.files(results_dir, pattern = "truths_cross_sectional")
all_truth_files <- all_truth_files[!grepl("all", all_truth_files)]
all_truth_lst <- lapply(as.list(all_truth_files), function(file) {
  readRDS(paste0(results_dir, file))
})
all_truths <- data.table::rbindlist(all_truth_lst, fill = TRUE) %>% 
  mutate(designation = gsub("-auc", "-autc", designation)) %>% 
  mutate(corr_between = ifelse(is.na(corr_between), 0, corr_between),
         corr_within = ifelse(is.na(corr_within), 0.5, corr_within))

all_output2 <- all_output_2 %>% 
  left_join(all_truths, by = c("corr_within", "corr_between", "designation", "varset"))
saveRDS(na.omit(all_output2, cols = 1), file = paste0(results_dir, "all_output.rds"))
saveRDS(na.omit(all_output2, cols = 1), file = paste0(here::here("..", "results", "sims"), "all_output.rds"))
