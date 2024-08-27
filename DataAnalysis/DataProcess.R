library(tidyverse)
library(ggplot2)
library(ggpubr)

# Import the raw data
# data <- read.table(here::here("DataAnalysis","Data","full2dat.txt"), header = TRUE)
# 
# # Viral load detection limit is 2
# # Define the viral rebound as two consecutive viral rise
# data <- data %>% 
#   mutate(`CD4` = ifelse(cd4 > 200, 1, 0)) %>% 
#   mutate(censor = ifelse(censor == -1, 1, 0)) %>% 
#   mutate(limit =0) %>% 
#   group_by(patid) %>% 
#   mutate(prev = lag(rna)) %>% 
#   mutate(diff = rna - prev) %>% 
#   mutate(nextdiff = lead(diff)) %>% 
#   mutate(rebound = ifelse(nextdiff > 0 & diff > 0, 1, 0)) %>% 
#   mutate(rebound = ifelse(is.na(rebound), 0, rebound)) %>% 
#   mutate(reboundtime = ifelse(rebound == 1, day, 0)) %>%
#   mutate(reboundtime = ifelse(sum(reboundtime) != 0, max(reboundtime), max(day)+1)) %>% 
#   mutate(reboundnew = ifelse(day >= reboundtime, 1, 0)) %>% 
#   dplyr::select(1:7, 13) %>% 
#   mutate(decay = ifelse(reboundnew == 1, 0, 1))
# 
# # Manually define some decay and rebound
# data <- data %>% 
#   mutate(decay = ifelse(patid == 610470, 1, decay)) %>% 
#   mutate(reboundnew = ifelse(decay == 1, 0, 1))
# data2 <- data
# data2$reboundnew[c(31, 53, 67, 68, 76, 112, 131, 148, 172, 189, 217, 314, 361)] <- 1
# data2 <- data2 %>% 
#   mutate(decay = ifelse(reboundnew == 1, 0, 1))
# 
# # Pivot long
# data.long <- pivot_longer(data2, cols=c(rna, CD4), names_to = "name", values_to = "observation")
# data.long <- data.long %>% 
#   mutate(observation = ifelse(is.na(observation), ".", observation)) %>% 
#   mutate(censor = ifelse(name == "CD4", 0, censor))
# 
# # New decay/rebound definition
# data3 <- data.long %>% group_by(patid) %>% mutate(reboundnew = ifelse(sum(reboundnew)!=0 & decay==1, lead(reboundnew, 2), reboundnew))
# data3 <- data3 %>% mutate(decay = ifelse(reboundnew == 1, 0, 1))
# 
# # Decay data
# decay <- data3 %>% 
#   filter(name == "rna") %>% 
#   filter(decay == 1) %>% 
#   mutate(observation = ifelse(censor == 1, log10(50), observation)) %>% 
#   dplyr::select(-c(cd4, limit, reboundnew, decay, name, censor))
# write.csv(decay, here::here("DataAnalysis","Data","Cleaned_data","full2dat_decay.csv"), row.names = FALSE)
# 
# # Rebound data
# rebound <- data3 %>% 
#   group_by(patid) %>% 
#   mutate(Ti = case_when(decay == 1 ~ day)) %>% 
#   mutate(Ti = max(Ti, na.rm = TRUE)) %>% 
#   filter(name == "rna") %>% 
#   filter(decay == 0) %>% 
#   mutate(observation = ifelse(censor == 1, log10(50), observation)) %>% 
#   mutate(time = day - Ti) %>% 
#   dplyr::select(-c(cd4, limit, reboundnew, decay, name, censor, day, Ti))
# write.csv(rebound, here::here("DataAnalysis","Data","Cleaned_data","full2dat_rebound.csv"), row.names = FALSE)
# 
# 
# # GLMM data
# glmm <- data3 %>% 
#   mutate(observation = ifelse(censor == 1, log10(50), observation)) %>% 
#   select(-c(limit, censor, reboundnew)) %>% 
#   pivot_wider(names_from = name, values_from = observation) %>% 
#   select(-c(cd4, decay))
# viral_mean <- mean(glmm$rna)
# viral_sd <- sd(glmm$rna)
# glmm <- glmm %>% 
#   mutate(viral_st = (rna - viral_mean)/viral_sd) %>% 
#   select(-rna)
# write.csv(rebound, here::here("DataAnalysis","Data","Cleaned_data","full2dat_glmm.csv"), row.names = FALSE)

filename = "raw_data.csv"
outliers = c(4,20,27)

######################### Clean raw data ================
# Read raw data
raw_data<-read.csv(here::here("DataAnalysis","Data","raw_data.csv"),header = T)
# Set the detection limit to be 40.
raw_data <- raw_data %>%
  mutate(DL = 40) 
# Find lowest detection limit for each patient
# raw_data <- raw_data %>%
#   group_by(PATIENT) %>%
#   mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
# If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
# Take log10 to RNA_V. Convert time from days to months. 
# Set lower boundary (limit) for the censoring interval to 0.
# Remove rows with unobserved RNA_V and time. 
raw_data <- raw_data %>% 
  mutate(censor = ifelse(RNA_V<=40, 1, 0)) %>% 
  mutate(RNA_V = ifelse(censor==1, 1/2*DL, RNA_V)) %>% 
  mutate(log10 = log10(RNA_V)) %>% 
  mutate(time_raw = days_from_seroco/30) %>% 
  mutate(limit = 0) %>% 
  filter(!is.na(RNA_V) & !is.na(time_raw))
# Remove outliers if there is any
if(!is.null(outliers)){
  raw_data <- raw_data %>% 
    filter(!(PATIENT %in% outliers))
}
# Find interruption time Ti for each patient. 
# Create a new variable "phase" to indicate viral decay and viral rebound.
# Change time for rebound data to the time since ART interruption 
# For example: if during ART time is c(1,2,3), viral rebound time is c(3.5,4,5), then 
#             the time since ART interruption for after_ART is c(0.5,1,2).
# Select useful variables.
data <- raw_data %>% 
  group_by(PATIENT) %>% 
  mutate(Ti = case_when(treatment == 1 ~ time_raw)) %>%  
  mutate(Ti = max(Ti, na.rm=TRUE)) %>% 
  mutate(phase = ifelse(time_raw <= Ti, "decay", "rebound")) %>% 
  mutate(time = ifelse(phase=="decay", time_raw, time_raw-Ti)) %>% 
  mutate(decay = ifelse(phase == "decay", 1, 0)) %>% 
  mutate(rebound = ifelse(phase == "rebound", 1, 0))
# Remove outliers
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
  filter(!(PATIENT == 32 & days_from_seroco ==  810)) #%>% 
  # filter(!(PATIENT == 13 & decay == 0)) %>% 
  # filter(!(PATIENT == 61 & decay == 0)) %>% 
  # filter(!(PATIENT == 73 & decay == 0)) %>% 
  # filter(!(PATIENT == 5 & decay == 0)) %>% 
  # filter(!(PATIENT == 14 & decay == 0)) %>% 
  # filter(!(PATIENT == 26 & decay == 0)) %>% 
  # filter(!(PATIENT == 6 & decay == 0)) %>% 
  # filter(!(PATIENT == 37 & decay == 0)) %>% 
  # filter(!(PATIENT == 41 & decay == 0)) 

# Split before_ART and after_ART based on phase.
before_ART <- data %>% filter(phase=="decay") %>% 
  dplyr::select(PATIENT, time, log10)
before_ART2 <- data %>% filter(phase=="decay") %>% 
  dplyr::select(PATIENT, time, log10, CD4_V) %>% 
  filter(!is.na(CD4_V))
baseline_CD4 <- aggregate(before_ART2$CD4_V, by=list(before_ART2$PATIENT), FUN=first)
colnames(baseline_CD4) <- c("PATIENT", "CD4")
before_ART2 <- merge(before_ART2, baseline_CD4, by = "PATIENT") %>% select(-CD4_V)
CD4.mean <- mean(baseline_CD4$CD4)
CD4.sd <- sd(baseline_CD4$CD4)
before_ART2 <- before_ART2 %>% 
  mutate(CD4 = (CD4-CD4.mean)/CD4.sd)

after_ART <- data %>% filter(phase=="rebound")%>% 
  dplyr::select(PATIENT, time, log10)

GLMM <- data %>% 
  mutate(`CD4` = ifelse(CD4_V > 500, 1, 0)) %>% 
  dplyr::select(PATIENT, days_from_seroco, log10, CD4) %>% 
  filter(!is.na(CD4)) 
viral_mean <- mean(GLMM$log10)
viral_sd <- sd(GLMM$log10)
GLMM <- GLMM %>% 
  mutate(viral_st = (log10-viral_mean)/viral_sd) %>% 
  mutate(time = days_from_seroco/30) %>% 
  select(-log10, -days_from_seroco)
GLMM$CD4 <- as.factor(GLMM$CD4)


GLMM2 <- data %>% 
  dplyr::select(PATIENT, days_from_seroco, censor, CD4_V) %>% 
  filter(!is.na(CD4_V)) 
CD4.mean <- mean(GLMM2$CD4_V)
CD4.sd <- sd(GLMM2$CD4_V)
GLMM2 <- GLMM2 %>% 
  mutate(CD4_st = (CD4_V-CD4.mean)/CD4.sd) %>% 
  mutate(time = days_from_seroco/30) %>% 
  select(-CD4_V, -days_from_seroco)

GLMM3 <- data %>% 
  dplyr::select(PATIENT, days_from_seroco, censor, CD4_V) %>% 
  filter(!is.na(CD4_V)) 
CD4.mean <- mean(GLMM3$CD4_V)
CD4.sd <- sd(GLMM3$CD4_V)
GLMM3 <- GLMM3 %>% 
  mutate(CD4_st = (CD4_V-CD4.mean)/CD4.sd) %>% 
  mutate(time = days_from_seroco/30) %>% 
  select(-CD4_V, -days_from_seroco) %>% 
  mutate(time2 = time^2)


# Save the two datasets locally. This is required when run Monolix.
write.csv(before_ART,here::here("DataAnalysis","Data","Cleaned_data","decay.csv"),row.names = FALSE)
write.csv(before_ART2,here::here("DataAnalysis","Data","Cleaned_data","decay2.csv"),row.names = FALSE)
write.csv(after_ART,here::here("DataAnalysis","Data","Cleaned_data","rebound.csv"),row.names = FALSE)
write.csv(GLMM,here::here("DataAnalysis","Data","Cleaned_data","GLMM.csv"),row.names = FALSE)
write.csv(GLMM2,here::here("DataAnalysis","Data","Cleaned_data","GLMM2.csv"),row.names = FALSE)
write.csv(GLMM3,here::here("DataAnalysis","Data","Cleaned_data","GLMM3.csv"),row.names = FALSE)
