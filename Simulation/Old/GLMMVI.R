library(lixoftConnectors)
library(mvtnorm)
library(nlme)
library(saemix)
library(tidyverse)
library(lme4)
library(MASS)
library(xtable)
library(berryFunctions)
library(kableExtra)
initializeLixoftConnectors(software="monolix")


##=====================Setting======================

num_sim <- 500  #number of simulation
num_patient <- 50  #number of patients

#Choose time from true_value dataset. 
# time <- c(0.5, 1.67, 3, 4.6, 6) #PATIENT #1
time <- c(0.5, 1.67, 3, 4.6, 6, 7.8, 9, 10.3, 11.2, 12.6)

#set true values for alpha0, alpha1
alpha0 <- 8
alpha1 <- -2
alpha2 <- 2
G <- diag(c(2, 0.5, 0.5))

#Create a data frame for simulation estimates.
estimates_SAEM <- estimates_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(estimates_SAEM) <- colnames(estimates_LME4) <- c("alpha0", "alpha1", "alpha2")

SE_SAEM <- SE_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(SE_SAEM) <- colnames(SE_LME4) <- c("alpha0", "alpha1", "alpha2")

CI_SAEM <- CI_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(CI_SAEM) <- colnames(CI_LME4) <- c("alpha0", "alpha1", "alpha2")

time_SAEM <- time_LME4 <- c(rep(NA, num_sim))

i=1

for (i in 1:num_sim){
  print(glue::glue("i=", i))
  ##==================Viral load before ART interruption======================
  
  #Create dataset for viral load before ART interruption: 1/2*num_patient patients with time1 and 1/2*num_patient patients with time2
  PATIENT <- rep(1:num_patient, each = length(time))
  day <- rep(time, num_patient)
  data <- cbind(PATIENT, day)
  
  # simulate b_i and e_ij
  # simulate a_i
  ai_sim <- mvrnorm(num_patient, c(0,0,0), G)
  colnames(ai_sim)=c("a1_sim","a2_sim","a3_sim")
  ai_sim <- cbind(PATIENT=c(1:num_patient),ai_sim)
  data <- merge(data, ai_sim, by="PATIENT")

  data_sim <- data %>% 
    mutate(x = rnorm(nrow(data), 0, 1)) %>% 
    mutate(xb = (alpha0 + a1_sim) + (alpha1 + a2_sim) * day + (alpha2 + a3_sim) * x) %>% 
    mutate(p = 1/(1+exp(-xb))) %>% 
    mutate(y = rbinom(n = nrow(data), size = 1, p = p)) %>% 
    dplyr::select(-a1_sim, -a2_sim, -a3_sim, -xb, -p)
  summary(as.factor(data_sim$y))

  # save the simulated data locally
  write.table(data_sim,"Data/data_GLMMVI.txt", sep = "," ,quote = FALSE,row.names = FALSE)
  
  # Fit viral decay model using Monolix.
  start.time <- Sys.time()
  data = list(dataFile = paste0('Data/data_GLMMVI.txt'),
              headerTypes =c("id","time", "regressor","observation"),
              observationTypes = list(y ="categorical"))
  modelFile = paste0('Model/model_GLMM3.txt')
  newProject(data = data, modelFile =modelFile)
  # getObservationInformation()
  # setErrorModel(viral = "constant")
  # setObservationDistribution(viral = "normal")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE)
  setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal")
  setPopulationParameterInformation(alpha0_pop = list(initialValue = 8), 
                                    alpha1_pop = list(initialValue = -2),
                                    alpha2_pop = list(initialValue = 2))

  setScenario(scenario)
  # run the estimation
  runScenario()
  
  # store the estimates in table
  estimates_SAEM[i,] <- getEstimatedPopulationParameters()[1:3]
  SE_SAEM[i,] <- getEstimatedStandardErrors()$stochasticApproximation[1:3,2]
  end.time <- Sys.time()
  time_SAEM[i] <- end.time - start.time
 
  
  start.time <- Sys.time()
  tryCatch({
    # cat("Iteration:", i, "\n")
    data_lme4 <- groupedData(y ~ day | PATIENT, data = data_sim)
    # Fit the generalized linear mixed model using glmer
    model_lme4 <- glmer(
      y ~ day + x + (1 + day + x|PATIENT),
      data = data_lme4, family = binomial, 
      start = list(fixef = c(8, -2, 2)))
    
    estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
    SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
  }, warning = function(w) {estimates_LME4[i, ] <<- SE_LME4[i, ] <<- NA}, 
  error = function(e) {estimates_LME4[i, ] <<- SE_LME4[i, ] <<-NA})
  end.time <- Sys.time()
  time_LME4[i] <- end.time - start.time
  
  write.csv(estimates_SAEM, "Results/GLMMVI/estimates_SAEM.csv", row.names = FALSE)
  write.csv(estimates_LME4, "Results/GLMMVI/estimates_LME4.csv", row.names = FALSE)
  write.csv(SE_SAEM, "Results/GLMMVI/SE_SAEM.csv", row.names = FALSE)
  write.csv(SE_LME4, "Results/GLMMVI/SE_LME4.csv", row.names = FALSE)
  write.csv(time_SAEM, "Results/GLMMVI/time_SAEM.csv", row.names = FALSE)
  write.csv(time_LME4, "Results/GLMMVI/time_LME4.csv", row.names = FALSE)
}  




#======================Calculate MSE and bias===================
# Read the results file
estimates_SAEM <- read.csv("Results/GLMMVI/estimates_SAEM.csv")
estimates_LME4 <- read.csv("Results/GLMMVI/estimates_LME4.csv")
SE_SAEM <- read.csv("Results/GLMMVI/SE_SAEM.csv")
SE_LME4 <- read.csv("Results/GLMMVI/SE_LME4.csv")
time_SAEM <- read.csv("Results/GLMMVI/time_SAEM.csv")$x
time_LME4 <- read.csv("Results/GLMMVI/time_LME4.csv")$x

true_value <- c(alpha0, alpha1, alpha2)

# Results table for SAEM 
# is.na(SE_SAEM) <- SE_SAEM == "NaN"
# SE_SAEM <- na.omit(SE_SAEM)
# SE_SAEM <- SE_SAEM %>% mutate_if(is.character,as.numeric)
# rownumber <- as.numeric(rownames(SE_SAEM))
# estimates_SAEM <- estimates_SAEM[rownumber, ]
# time_SAEM <- time_SAEM[rownumber]
for (i in 1:num_sim){
  CI_SAEM[i, ] = (estimates_SAEM[i, ] - 1.96 * SE_SAEM[i, ] <= true_value & estimates_SAEM[i, ] + 1.96 * SE_SAEM[i, ] >= true_value)
}
CI_SAEM <- na.omit(CI_SAEM)
rMSE_SAEM <- bias_SAEM <- bias_perc_SAEM <- vector("double", length(true_value))
for (i in seq_along(true_value)){
  rMSE_SAEM[i]=sqrt(sum((estimates_SAEM[,i]-true_value[i])^2)/nrow(estimates_SAEM))/true_value[i]
  bias_SAEM[i]=sum(estimates_SAEM[,i]-true_value[i])/nrow(estimates_SAEM)
  bias_perc_SAEM[i]=bias_SAEM[i]/true_value[i]
}
Coverage_SAEM <- sapply(CI_SAEM, mean)
result_SAEM <- cbind(true_value, estimates = sapply(estimates_SAEM, mean), se_m = sapply(SE_SAEM, mean), se_s = sapply(estimates_SAEM, sd), rMSE_SAEM,bias_SAEM,bias_perc_SAEM, Coverage_SAEM, NC = c(num_sim-nrow(estimates_SAEM), NA, NA), time = c(sum(time_SAEM), NA, NA))
rownames(result_SAEM)=c("alpha0", "alpha1", "alpha2")

# Results table for LME4 
estimates_LME4 <- na.omit(estimates_LME4)
SE_LME4 <- na.omit(SE_LME4)
time_LME4 <- na.omit(time_LME4)
for (i in 1:nrow(estimates_LME4)){
  CI_LME4[i, ] = (estimates_LME4[i, ] - 1.96 * SE_LME4[i, ] <= true_value & estimates_LME4[i, ] + 1.96 * SE_LME4[i, ] >= true_value)
}
CI_LME4 <- na.omit(CI_LME4)
rMSE_LME4 <- bias_LME4 <- bias_perc_LME4 <- vector("double", length(true_value))
for (i in seq_along(true_value)){
  rMSE_LME4[i]=sqrt(sum((estimates_LME4[,i]-true_value[i])^2)/nrow(estimates_LME4))/true_value[i]
  bias_LME4[i]=sum(estimates_LME4[,i]-true_value[i])/nrow(estimates_LME4)
  bias_perc_LME4[i]=bias_LME4[i]/true_value[i]
}
Coverage_LME4 <- sapply(CI_LME4, mean)
result_LME4 <- cbind(true_value,estimates = sapply(estimates_LME4, mean), se_m = sapply(SE_LME4, mean), se_s = sapply(estimates_LME4, sd), rMSE_LME4,bias_LME4,bias_perc_LME4, Coverage_LME4, NC = c(num_sim - nrow(estimates_LME4), NA, NA), time = c(sum(time_LME4), NA, NA))
rownames(result_LME4)=c("alpha0", "alpha1", "alpha2")

# Latex Table
table_names <- c( "Method", "Time (s)", "NC","Parameter","True Value", "Estimate", "SE_M",  "SE_S", "Bias (%)","rMSE (%)","Coverage")
table <- data.frame(matrix(nrow=6, ncol=length(table_names)))
colnames(table) <- table_names
table$Parameter <- c("$\\alpha_0$","$\\alpha_1$","$\\alpha_2$","$\\alpha_0$","$\\alpha_1$","$\\alpha_2$")
table$`True Value` <- c(8, -2, 2, 8, -2, 2)
table$Method <- c("SAEM",NA,NA,"LME4",NA,NA)
table$Estimate <- c(result_SAEM[, 2], result_LME4[,2])
table$SE_M <- c(result_SAEM[, 3],  result_LME4[,3])
table$SE_S <- c(result_SAEM[, 4], result_LME4[,4])
table$`Bias (%)` <- c(result_SAEM[, 7]*100, result_LME4[,7]*100)
table$`rMSE (%)` <- c(result_SAEM[, 5]*100, result_LME4[,5]*100)
table$Coverage <-  c(result_SAEM[, 8]*100, result_LME4[,8]*100)
table$`Time (s)` <- c(result_SAEM[, 10], result_LME4[,10])
table$NC <- c(result_SAEM[, 9], result_LME4[,9])

latex_table<-xtable(table, type = "latex",align=c("cccccccccccc"))
digits(latex_table)<-c(0,2,1,0,1,1,2,2,2,2,2,2)
print(latex_table, file = "Table/GLMMVI.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,6))



