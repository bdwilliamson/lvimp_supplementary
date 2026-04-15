# summarize cohort study analysis dataset for IMATS VIM analyses based on SRS3
library("here")
library("dplyr")
library("gt")
library("gtsummary")
library("bstfun")
library("kableExtra")
library("flextable")
library("huxtable")
library("tidyr")

# read in the dataset ----------------------------------------------------------
data_dir <- "<the directory where the data are stored>"
results_dir <- "<the directory where the results are stored>"
analysis_dataset <- readRDS(paste0(data_dir, "imats_srs3_cohort_study_analysis_dataset.rds"))

# create tables with summary statistics ----------------------------------------
nrow(analysis_dataset)
ncol(analysis_dataset)
X <- analysis_dataset %>% 
  select(-study_id, -visit_n, -timepoint, -event90)
ncol(X)

# number of unique people overall; number of people with each number of visit
length(unique(analysis_dataset$study_id))
analysis_dataset %>%
  group_by(study_id) %>%
  summarize(n_visits = n()) %>%
  ungroup() %>%
  pull(n_visits) %>%
  table()

# back-engineer categorical variables for nicer display
table_1_data <- analysis_dataset %>% 
  mutate(site_fct = factor(case_when(site_hpi == 1 ~ "HealthPartners",
                                     site_hpi == 0 ~ "KPWA"),
                           levels = c("HealthPartners", "KPWA"),
                           labels = c("HealthPartners", "KPWA")),
        item9_fct = factor(item9_visit,
                            levels = c(-5, 0, 1, 2, 3),
                            labels = c("Not measured", "0", "1", "2", "3")),
        phq8_fct = factor(cut(phq8_visit, breaks = c(-10, -1, 4, 10, 15, 20, 25), labels = c("Not measured", "0--4", "5--10", "11--15", "16--20", "21 or higher"))),
         age_fct = factor(cut(age, breaks = c(10, 17, 29, 44, 64, 200), labels = c("11--17", "18--29", "30--44", "45--64", "65 and older"))),
         sex_fct = factor(case_when(sex_f == 1 ~ "Female",
                                    sex_u == 1 ~ "Unknown sex",
                                    sex_f == 0 & sex_o == 0 & sex_u == 0 ~ "Male"),
                          levels = c("Female", "Male", "Unknown sex"),
                          labels = c("Female", "Male", "Unknown sex")),
         race_fct = factor(case_when(race_as == 1 ~ "Asian",
                                     race_ba == 1 ~ "Black or African American",
                                     race_hp == 1 ~ "Native Hawaiian or Pacific Islander",
                                     race_in == 1 ~ "American Indian/Alaska Native",
                                     race_mu == 1 ~ "Multiple races",
                                     race_ot == 1 ~ "Other race",
                                     race_un == 1 ~ "Unknown race",
                                     race_as == 0 & race_ba == 0 & race_hp == 0 & race_in == 0 & 
                                       race_mu == 0 & race_ot == 0 & race_un == 0 ~ "White"),
                           levels = c("Asian", "Black or African American", "Native Hawaiian or Pacific Islander",
                                      "American Indian/Alaska Native", "Multiple races", "Other race", "Unknown race",
                                      "White"),
                           labels = c("Asian", "Black or African American", "Native Hawaiian or Pacific Islander",
                                      "American Indian/Alaska Native", "Multiple races", "Other race", "Unknown race",
                                      "White")),
         income_fct = factor(case_when(income_unknown == 1 ~ "Unknown income",
                                       income_lt25 == 1 ~ "Income < $25K",
                                       income_lt40 == 1 ~ "$25K <= Income < $40K",
                                       income_unknown == 0 & income_lt25 == 0 & income_lt40 == 0 ~ "Income >= $40K"),
                             levels = c("Unknown income", "Income < $25K", "$25K <= Income < $40K", "Income >= $40K"),
                             labels = c("Unknown income", "Income < $25K", "$25K <= Income < $40K", "Income >= $40K")),
         edu_fct = factor(case_when(coll_deg_unknown == 1 ~ "Percentage unknown",
                                    coll_deg_lt25 == 1 ~ "<25%",
                                    coll_deg_unknown == 0 & coll_deg_lt25 == 0 ~ ">=25%"),
                          levels = c("Percentage unknown", "<25%", ">=25%"),
                          labels = c("Percentage unknown", "<25%", ">=25%")),
         ins_fct = factor(case_when(ins_aca_num == 2 ~ "Affordable Care Act",
                                    ins_medicaid_num == 2 ~ "Medicaid",
                                    ins_medicare_num == 2 ~ "Medicare",
                                    ins_commercial_num == 2 ~ "Commercial",
                                    ins_privatepay_num == 2 ~ "Private pay",
                                    ins_statesubsidized_num == 2 ~ "State-subsidized",
                                    ins_selffunded_num == 2 ~ "Self-funded",
                                    ins_highdeductible_num == 2 ~ "High-deductible",
                                    ins_other_num == 2 ~ "Other coverage",
                                    ins_aca_num == 0 & ins_medicaid_num == 0 & ins_medicare_num == 0 & 
                                      ins_commercial_num == 0 & ins_privatepay_num == 0 & ins_statesubsidized_num == 0 & 
                                      ins_selffunded_num == 0 & ins_highdeductible_num == 0 & ins_other_num == 0 ~ "Unknown coverage",
                                    ins_aca_num == 1 & ins_medicaid_num == 1 & ins_medicare_num == 1 & 
                                      ins_commercial_num == 1 & ins_privatepay_num == 1 & ins_statesubsidized_num == 1 & 
                                      ins_selffunded_num == 1 & ins_highdeductible_num == 1 & ins_other_num == 1 ~ "No coverage"),
                          levels = c("Affordable Care Act", "Medicaid", "Medicare", "Commercial", "Private pay", "State-subsidized",
                                     "Self-funded", "High-deductible", "Other coverage", "Unknown coverage", "No coverage"),
                          labels = c("Affordable Care Act", "Medicaid", "Medicare", "Commercial", "Private pay", "State-subsidized",
                                     "Self-funded", "High-deductible", "Other coverage", "Unknown coverage", "No coverage")),
         visit_type_fct = factor(case_when(visit_pc1 == 1 ~ "Primary care (SRS2 diagnosis categories)",
                                           visit_pc2 == 1 ~ "Primary care (additional categories)",
                                           visit_pc1 == 0 & visit_pc2 == 0 ~ "Mental health specialty"),
                                 levels = c("Primary care (SRS2 diagnosis categories)", "Primary care (additional categories)", "Mental health specialty"),
                                 labels = c("Primary care (SRS2 diagnosis categories)", "Primary care (additional categories)", "Mental health specialty"))
  ) %>%
  mutate(across(ends_with("tot_days_last12m"), .fns = ~ as.numeric(.x > 0), .names = "{.col}_any"))
  
# visit-process summaries; this is only for internal use
table_1_visit_process <- table_1_data %>% 
  select(timepoint, visit_n, days_since_visit1, mths_since_prev, enr_calc,
         visit_type_fct) %>% 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      visit_n ~ "Visit sequence number",
      days_since_visit1 ~ "Days since first qualifying visit",
      mths_since_prev ~ "Months since previous visit",
      enr_calc ~ "Months of continuous pre-visit enrollment (trunc. at 60)",
      visit_type_fct ~ "Visit type"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA,
                             stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA)
# demographic variables
table_1_demographics <- table_1_data %>% 
  select(timepoint, age_fct, sex_fct) |> 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_fct ~ "Age in years",
      sex_fct ~ "Sex"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA)
# race and ethnicity variables
table_1_race <- table_1_data %>% 
  select(timepoint, hispanic, race_fct) %>% 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      hispanic ~ "Hispanic ethnicity",
      race_fct ~ "Self-reported race"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) %>%
  gtsummary::modify_table_styling(columns = label, rows = label == "Self-reported race",
                                  footnote = "Individuals who reported more than one listed race or ethnicity contribute to all selected racial and ethnic subgroups.")

# diagnosis and prescribing variables: presence or absence in the last 12 months
table_1_dx <- table_1_data %>% 
  select(timepoint, ends_with("tot_days_last12m_any")) %>% 
  select(-starts_with("asa"), -starts_with("lsa"), -starts_with("osa"), -starts_with("aip"), -starts_with("mh")) %>%
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      add_tot_days_last12m_any ~ "Atention deficit disorder",
      alc_tot_days_last12m_any ~ "Alcohol use disorder",
      anx_tot_days_last12m_any ~ "Anxiety",
      asd_tot_days_last12m_any ~ "Autism spectrum disorder",
      ast_tot_days_last12m_any ~ "Asthma",
      con_tot_days_last12m_any ~ "Conduct disorder",
      dem_tot_days_last12m_any ~ "Dementia/cognitive disorder",
      dep_tot_days_last12m_any ~ "Depression",
      dia_tot_days_last12m_any ~ "Diabetes",
      dru_tot_days_last12m_any ~ "Drug use disorder",
      eat_tot_days_last12m_any ~ "Eating disorder",
      opd_tot_days_last12m_any ~ "Psychosis",
      pai_tot_days_last12m_any ~ "Pain",
      per_tot_days_last12m_any ~ "Personality disorder",
      pts_tot_days_last12m_any ~ "Post-traumatic stress disorder",
      tbi_tot_days_last12m_any ~ "Traumatic brain injury"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) %>%
  gtsummary::modify_table_styling(columns = label, variable %in% paste0(c("add", "alc", "anx", "asd", "ast", "con", "dem", "dep", "dia", "dru", "eat", "opd", "pai", "per", "pts", "tbi"), "_tot_days_last12m_any"), 
                                  footnote = "At least one diagnosis in the past 12 months.")

table_1_prior_self_harm <- table_1_data %>%
  select(timepoint, ends_with("tot_days_last12m_any")) %>%
  select(timepoint, starts_with("asa")) |> 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      asa_tot_days_last12m_any ~ "Prior self-harm"
    )
  ) %>%
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) %>%
  gtsummary::modify_table_styling(columns = label, rows = variable %in% paste0(c("asa", "lsa", "osa", "aip"), "_tot_days_last12m_any"), 
                                  footnote = "At least one injury, poisoning, or attempt in the past 12 months.")

table_1_encounters <- table_1_data %>%
  select(timepoint, ends_with("tot_days_last12m_any")) %>%
  select(timepoint, starts_with("mh")) %>%
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      mhe_tot_days_last12m_any ~ "Mental health emergency / urgent care encounter",
      mhi_tot_days_last12m_any ~ "Mental health inpatient encounter",
      mho_tot_days_last12m_any ~ "Mental health outpatient encounter"
    )
  ) %>%
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) %>%
  gtsummary::modify_table_styling(columns = label, rows = variable %in% paste0(c("mhe", "mhi", "mho"), "_tot_days_last12m_any"), 
                                  footnote = "At least one encounter in the past 12 months.")

# PHQ and charlson variables, event
table_1_event <- table_1_data %>%
  select(timepoint, event90) %>%
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      event90 ~ "Suicide attempt"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) |> 
  gtsummary::modify_table_styling(columns = label, rows = label == "Suicide attempt", 
                                  footnote = "Suicide attempt in the 90 days following the mental health care visit.")
table_1_phq8 <- table_1_data |> 
  select(timepoint, phq8_fct) |> 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      phq8_fct ~ "PHQ-8 total score"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) 

table_1_phq9 <- table_1_data |> 
  select(timepoint, item9_fct) |> 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      item9_fct ~ "PHQi9"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA) |> 
  gtsummary::modify_table_styling(columns = label, rows = label == "PHQi9",
                                  footnote = "PHQi9 and PHQ-8 total score measured on the day of the mental health visit")
table_1_charlson <- table_1_data %>% 
  select(timepoint, charlson_score) %>% 
  gtsummary::tbl_summary(
    by = "timepoint",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      charlson_score ~ "Charlson comorbidity index"
    )
  ) %>% 
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Timepoint (days from initial visit)**") %>%
  gtsummary::add_overall(col_label = "**Overall**\nN = {N}") %>%
  gtsummary::modify_footnote(stat_0 ~ NA, stat_1 = NA, stat_2 = NA, stat_3 = NA, stat_4 = NA, stat_5 = NA, stat_6 = NA)

# final table
table_1_internal <- gtsummary::tbl_stack(
 list(table_1_event, table_1_demographics, table_1_race, table_1_phq9, table_1_phq8, table_1_charlson,
       table_1_dx, table_1_prior_self_harm, table_1_encounters, table_1_visit_process),
  group_header = c("Outcome", "Demographic\nvariables", "Race and ethnicity\nvariables", "", "", "",
                   "Diagnosis\nvariables", "Prior self-harm\nvariables", "Encounter\nvariables",
                   "Visit process\nvariables"),
  quiet = TRUE
)
flextable::save_as_docx(path = paste0(results_dir, "table_1_internal.docx"),
                        table_1_internal %>% gtsummary::as_flex_table())

table_1 <- gtsummary::tbl_stack(
   list(table_1_event, table_1_demographics, table_1_prior_self_harm, table_1_phq9, table_1_phq8),
  group_header = c("", "Demographic\n variables", "", "", "")
) %>% 
  as_gt() %>%
  opt_footnote_marks(marks = "letters") %>%
  tab_caption(caption = "Cohort description for sample used to estimate variable importance in suicide attempt risk prediction models\\label{tab:table_1}.") %>%
  tab_options(table.font.size = 2)
table_1_latex <- as.character(table_1 %>% as_latex())
# hack to get width to work
y <- unlist(strsplit(table_1_latex, "\n\\\\"))
y2 <- y
y2_header <- y2[5]
y2_header_split <- unlist(strsplit(y2_header, "&"))
y2_header_new <- gsub("\n", " \\\\newline ", y2_header_split)
y2[5] <- paste0(y2_header_new, collapse = "& ")
y2[2] <- "begin{longtable}{p{1in}p{0.70in}p{0.68in}p{0.68in}p{0.68in}p{0.68in}p{0.68in}p{0.68in}p{0.68in}}"
y3 <- c(y2[1:2], "caption{Cohort description for sample used to estimate variable importance in suicide attempt risk prediction models.}\\label{tab:table_1}\\\\", y2[3:length(y2)][-c(4:5)]) # -4,5 gets rid of a "midrule" and an empty row
table_1_latex2 <- paste0(y3, collapse = "\n\\")
write(table_1_latex2, paste0(results_dir, "table_1.tex"))
