library(tidyverse)
library(ggplot2)
library(ggpubr)

filename = "raw_data.csv"
outliers = c(4,20,27)

######################### Clean raw data ================
# Read raw data
raw_data <- read.csv(here::here("DataAnalysis","Data","raw_data.csv"),header = T)
# Set the detection limit to be 40.
raw_data <- raw_data %>%
  mutate(DL = 40) 
# Data processing
raw_data <- raw_data %>% 
  mutate(censor = ifelse(RNA_V<=40, 1, 0)) %>% # Create a binary censoring indicator variable `censor`
  mutate(RNA_V = ifelse(censor==1, 1/2*DL, RNA_V)) %>% # If `censor` equals 1, replace `RNA_V` with half DL
  mutate(log10 = log10(RNA_V)) %>% # Compute the log10 transformation of `RNA_V`
  mutate(time_raw = days_from_seroco/30) %>% # Convert `days_from_seroco` to months by dividing by 30
  mutate(limit = 0) %>% # Create a new variable `limit` and set its value to 0 for all observations.
  filter(!is.na(RNA_V) & !is.na(time_raw)) # Filter the dataset to remove any rows where either `RNA_V` or `time_raw` is missing
# Remove outliers if there is any
if(!is.null(outliers)){
  raw_data <- raw_data %>% 
    filter(!(PATIENT %in% outliers))
}
# Data processing
data <- raw_data %>% 
  group_by(PATIENT) %>% 
  mutate(Ti = case_when(treatment == 1 ~ time_raw)) %>% # If `treatment` equals 1, assign `time_raw` to `Ti
  mutate(Ti = max(Ti, na.rm=TRUE)) %>% # Replace `Ti` with the maximum value of `Ti` within each patient group
  mutate(phase = ifelse(time_raw <= Ti, "decay", "rebound")) %>% # If `time_raw` is less than or equal to `Ti`, assign "decay" to `phase`; otherwise, assign "rebound"
  mutate(time = ifelse(phase=="decay", time_raw, time_raw-Ti)) %>% # If in the "decay" phase, keep `time_raw`; if in the "rebound" phase, subtract `Ti` from `time_raw` to reset the time counter from the start of the rebound phase.
  mutate(decay = ifelse(phase == "decay", 1, 0)) %>% # Create a binary indicator `decay`
  mutate(rebound = ifelse(phase == "rebound", 1, 0)) # Create a binary indicator `rebound`
# Remove outliers manually
data <- data %>% 
  filter(!(PATIENT == 4 & phase == "rebound")) %>% 
  filter(!(PATIENT == 30 & phase == "rebound")) %>% 
  filter(!(PATIENT == 5 & days_from_seroco == 1089)) %>% 
  filter(!(PATIENT == 11 & days_from_seroco == 712)) %>% 
  filter(!(PATIENT == 11 & days_from_seroco == 748)) %>% 
  filter(!(PATIENT == 12 & days_from_seroco == 989)) %>% 
  filter(!(PATIENT == 17 & days_from_seroco == 740)) %>% 
  filter(!(PATIENT == 27)) %>% 
  filter(!(PATIENT == 35 & days_from_seroco == 1181)) %>% 
  filter(!(PATIENT == 42 & days_from_seroco == 1190)) %>% 
  filter(!(PATIENT == 32 & days_from_seroco ==  742)) %>% 
  filter(!(PATIENT == 32 & days_from_seroco ==  810)) 

# Find the baseline CD4 value.
before_ART <- data %>% filter(phase=="decay") %>% 
  dplyr::select(PATIENT, time, log10, CD4_V) %>% 
  filter(!is.na(CD4_V))
baseline_CD4 <- aggregate(before_ART$CD4_V, by=list(before_ART$PATIENT), FUN=first)
colnames(baseline_CD4) <- c("PATIENT", "CD4")
before_ART <- merge(before_ART, baseline_CD4, by = "PATIENT") %>% dplyr::select(-CD4_V)

# Standardize the baseline CD4 value.
CD4.mean <- mean(baseline_CD4$CD4)
CD4.sd <- sd(baseline_CD4$CD4)
before_ART <- before_ART %>% 
  mutate(CD4 = (CD4-CD4.mean)/CD4.sd) 

# GLMM: whether viral load being below the detection limit is associated with CD4 values.
GLMM <- data %>% 
  dplyr::select(PATIENT, time_raw, censor, CD4_V) %>% 
  filter(!is.na(CD4_V)) 
CD4.mean <- mean(GLMM$CD4_V)
CD4.sd <- sd(GLMM$CD4_V)
GLMM <- GLMM %>% 
  mutate(CD4_st = (CD4_V-CD4.mean)/CD4.sd) %>% 
  dplyr::select(-CD4_V)

# Save the two datasets locally. 
write.csv(before_ART,here::here("DataAnalysis","Data","Cleaned_data","decay.csv"), row.names = FALSE)
write.csv(GLMM,here::here("DataAnalysis","Data","Cleaned_data","GLMM.csv"), row.names = FALSE)
