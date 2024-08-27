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
num_patient <- 200  #number of patients
# censor_value <- log10(40)  #censor value

#Set true values for p1, b1, p2, sigma1, and covariance matrix D for viral loads before ART interruption
#Follow model 1: y_ij=log10(p_1i+p_2i*exp(-lambda_1i*t_ij))+e_ij
#Random effects: p_1i=p1+b_1i,p_2i=p1+b_2i,lambda_1i=lambda1+b_3i,
#where e_ij~N(0,sigma1^2) and b_i~N(0,D)
sigma <- 0.1
p1 <- 17
lambda1 <- 4
p2 <- 3
D <- diag(c(2, 0.1, 0.1))

#Choose time from true_value dataset. 
# time <- c(0.5, 1.2, 1.67, 2.4, 3, 3.3, 4.2, 4.6, 5.3, 6) #PATIENT #1
time <- c(0.5, 0.8, 1.2, 1.67, 2.2, 2.4, 2.7, 3, 3.3, 3.6, 3.9, 4.2, 4.6, 5.3, 6) #PATIENT #1

#Create a data frame for simulation estimates.
estimates_SAEM <- estimates_NLME <- estimates_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(estimates_SAEM) <- colnames(estimates_NLME) <- colnames(estimates_LME4) <- c("p1","lambda1","p2")

SE_SAEM <- SE_NLME <- SE_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(SE_SAEM) <- colnames(SE_NLME) <- colnames(SE_LME4) <- c("p1","lambda1","p2")

CI_SAEM <- CI_NLME <- CI_LME4 <- data.frame(matrix(nrow=num_sim,ncol = 3))
colnames(CI_SAEM) <- colnames(CI_NLME) <- colnames(CI_LME4) <- c("p1","lambda1","p2")

time_SAEM <- time_NLME <- time_LME4 <- c(rep(NA, num_sim))

Model <- function(p1 ,b1 ,p2, t) {
  log10 ( exp(p1-b1*t)+ exp (p2))
}

ModelGradient <- deriv(
  body(Model)[[2]], 
  namevec = c("p1", "b1", "p2"), 
  function.arg = Model
)

i=1
fail = 0
for (i in 1:num_sim){
  print(glue::glue("i=", i))
  ##==================Viral load before ART interruption======================
  
  #Create dataset for viral load before ART interruption: 1/2*num_patient patients with time1 and 1/2*num_patient patients with time2
  PATIENT <- rep(1:num_patient, each = length(time))
  day <- rep(time, num_patient)
  data <- cbind(PATIENT, day)
  
  #simulate b_i and e_ij
  bi_sim <- mvrnorm(num_patient, c(0,0,0), D)
  colnames(bi_sim)=c("b1_sim","b2_sim","b3_sim")
  bi_sim <- cbind(PATIENT=c(1:num_patient),bi_sim)
  data <- merge(data, bi_sim, by="PATIENT")
  data <- data %>% 
    mutate(e = rnorm(nrow(data), 0, sigma))
  
  #simulate viral load and use censor_value to indicate censor. censor=1 if simulated viral load < censor_value.
  # data_saem <- data %>% 
  #   mutate(viral = log10(exp(p1+b1_sim-(lambda1+b3_sim)*day)+exp(p2+b2_sim))+e) %>%
  #   mutate(censor = ifelse(viral < censor_value, 1, 0)) %>% 
  #   mutate(viral = ifelse(viral < censor_value, censor_value, viral)) %>% 
  #   mutate(limit = 0) %>% 
  #   dplyr::select(-b1_sim, -b2_sim, -b3_sim, -e)
  # summary(data_saem$censor)
  data_sim <- data %>% 
    mutate(viral = log10(exp(p1+b1_sim-(lambda1+b3_sim)*day)+exp(p2+b2_sim))+e)%>% 
    dplyr::select(-b1_sim, -b2_sim, -b3_sim, -e)
  # # impute the censored value by half the detection limit for data_nlme
  # data_nlme <- data_saem %>% 
  #   mutate(viral = ifelse(censor == 1, 0.5*censor_value, viral)) %>% 
  #   dplyr::select(-censor, -limit)
  # ggplot(data_saem[PATIENT%in%c(1:5),], aes(x=day, y=viral)) +
  #   geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
  #   geom_line(aes(group=PATIENT)) +
  #   scale_x_continuous("Day") +
  #   scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  #   scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
  #   labs(color = "ART status", fill = "Data type")+
  #   # ggtitle("Plot for all observations")+
  #   theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  ggplot(data_sim[PATIENT%in%c(1:5),], aes(x=day, y=viral)) +
    geom_point() +
    geom_line(aes(group=PATIENT)) +
    scale_x_continuous("Time") +
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    # ggtitle("Plot for all observations")+
    theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  
  # save the simulated data locally
  write.table(data_sim,"Data/data_Decay_n50_ni15.txt", sep = "," ,quote = FALSE,row.names = FALSE)
  
  # Fit viral decay model using Monolix.
  start.time <- Sys.time()
  data = list(dataFile = paste0('Data/data_Decay_n50_ni15.txt'),
              headerTypes =c("id","time","observation"),
              observationTypes = list(viral="continuous"))
  modelFile = paste0('Model/model_decay.txt')
  newProject(data = data, modelFile =modelFile)
  # getObservationInformation()
  setErrorModel(viral = "constant")
  setObservationDistribution(viral = "normal")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE)
  setPopulationParameterInformation(p1_pop = list(initialValue = 17), 
                                    b1_pop = list(initialValue = 4),
                                    p2_pop = list(initialValue = 3))
  setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal")
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
    data_nlme <- groupedData(viral ~ day | PATIENT, data = data_sim)
    model_nlme <- nlme(viral ~ log10(exp(p1-b1*day)+exp(p2)),
                       fixed = p1+b1+p2 ~ 1,
                       random = pdDiag(list(p1 ~ 1, b1 ~ 1, p2 ~ 1)),
                       data = data_nlme,
                       start = c(p1 = 17, b1 = 4, p2 = 3))
    
    estimates_NLME[i, ] <- summary(model_nlme)$tTable[,1]
    SE_NLME[i, ] <- summary(model_nlme)$tTable[,2]
  }, warning = function(w) {estimates_NLME[i, ] <<- SE_NLME[i, ] <<- NA}, 
  message = function(m) {estimates_NLME[i, ] <<- SE_NLME[i, ] <<-NA}, 
  error = function(e) {estimates_NLME[i, ] <<- SE_NLME[i, ] <<- NA})
  
  end.time <- Sys.time()
  time_NLME[i] <- end.time - start.time
  
  start.time <- Sys.time()
  # Fit the non-linear mixed effects model
  # model_lme4 <- nlmer(
  #   # Response
  #   viral ~ 
  #     # Fixed effects
  #     ModelGradient(p1, b1, p2, day) ~ 
  #     # Random effects
  #     (p1 | PATIENT) + (b1 | PATIENT) + (p2 | PATIENT), 
  #   # Data
  #   data = data_sim, 
  #   start =  c(p1 = 17, b1 = 4, p2 = 3))
  # 
  # estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
  # SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
  # end.time <- Sys.time()
  # time_LME4[i] <- end.time - start.time
  
  nform <- ~log10 (exp(p1-b1*t)+ exp (p2))
  nfun <- deriv(nform, namevec = c("p1", "b1", "p2"), 
                function.arg = c("t", "p1", "b1", "p2"))
  tryCatch({
  model_lme4 <- nlmer(viral ~ nfun(day, p1, b1, p2) ~ (p1 | PATIENT) + (p2 | PATIENT) + (b1 | PATIENT), data_sim, 
                      start = c(p1 = 17, b1 = 4, p2 = 3))
  estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
  SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
  }, error = function(e) {estimates_LME4[i, ] <<- SE_LME4[i, ] <<- NA})
  end.time <- Sys.time()
  time_LME4[i] <- end.time - start.time
}  

write.csv(estimates_SAEM, "Results/Decay_n200_ni15/estimates_SAEM.csv", row.names = FALSE)
write.csv(estimates_NLME, "Results/Decay_n200_ni15/estimates_NLME.csv", row.names = FALSE)
write.csv(estimates_LME4, "Results/Decay_n200_ni15/estimates_LME4.csv", row.names = FALSE)
write.csv(SE_SAEM, "Results/Decay_n200_ni15/SE_SAEM.csv", row.names = FALSE)
write.csv(SE_NLME, "Results/Decay_n200_ni15/SE_NLME.csv", row.names = FALSE)
write.csv(SE_LME4, "Results/Decay_n200_ni15/SE_LME4.csv", row.names = FALSE)
write.csv(time_SAEM, "Results/Decay_n200_ni15/time_SAEM.csv", row.names = FALSE)
write.csv(time_NLME, "Results/Decay_n200_ni15/time_NLME.csv", row.names = FALSE)
write.csv(time_LME4, "Results/Decay_n200_ni15/time_LME4.csv", row.names = FALSE)


#======================Calculate MSE and bias===================
# Read the results file
estimates_SAEM <- read.csv("Results/Decay_n200_ni15/estimates_SAEM.csv")
estimates_NLME <- read.csv("Results/Decay_n200_ni15/estimates_NLME.csv")
estimates_LME4 <- read.csv("Results/Decay_n200_ni15/estimates_LME4.csv")
SE_SAEM <- read.csv("Results/Decay_n200_ni15/SE_SAEM.csv")
SE_NLME <- read.csv("Results/Decay_n200_ni15/SE_NLME.csv")
SE_LME4 <- read.csv("Results/Decay_n200_ni15/SE_LME4.csv")
time_SAEM <- read.csv("Results/Decay_n200_ni15/time_SAEM.csv")$x
time_NLME <- read.csv("Results/Decay_n200_ni15/time_NLME.csv")$x
time_LME4 <- read.csv("Results/Decay_n200_ni15/time_LME4.csv")$x

true_value <- c(p1, lambda1, p2)

# Results table for SAEM 
for (i in 1:num_sim){
  CI_SAEM[i, ] = (estimates_SAEM[i, ] - 1.96 * SE_SAEM[i, ] <= true_value & estimates_SAEM[i, ] + 1.96 * SE_SAEM[i, ] >= true_value)
}
rMSE_SAEM <- bias_SAEM <- bias_perc_SAEM <- vector("double", length(true_value))
for (i in seq_along(true_value)){
  rMSE_SAEM[i]=sqrt(sum((estimates_SAEM[,i]-true_value[i])^2)/num_sim)/true_value[i]
  bias_SAEM[i]=sum(estimates_SAEM[,i]-true_value[i])/num_sim
  bias_perc_SAEM[i]=bias_SAEM[i]/true_value[i]
}
Coverage_SAEM <- sapply(CI_SAEM, mean)
result_SAEM <- cbind(true_value, estimates = sapply(estimates_SAEM, mean), se_m = sapply(SE_SAEM, mean), se_s = sapply(estimates_SAEM, sd), rMSE_SAEM,bias_SAEM,bias_perc_SAEM, Coverage_SAEM, NC = c(0, NA,NA), time = c(sum(time_SAEM), NA, NA))
rownames(result_SAEM)=c("p1","lambda1","p2")

# Results table for NLME
for (i in 1:num_sim){
  CI_NLME[i, ] = (estimates_NLME[i, ] - 1.96 * SE_NLME[i, ] <= true_value & estimates_NLME[i, ] + 1.96 * SE_NLME[i, ] >= true_value)
}
CI_NLME <- na.omit(CI_NLME)
estimates_NLME <- na.omit(estimates_NLME)
SE_NLME <- na.omit(SE_NLME)
rMSE_NLME <- bias_NLME <- bias_perc_NLME <- vector("double", length(true_value))
for (i in seq_along(true_value)){
  rMSE_NLME[i]=sqrt(sum((estimates_NLME[,i]-true_value[i])^2)/nrow(estimates_NLME))/true_value[i]
  bias_NLME[i]=sum(estimates_NLME[,i]-true_value[i])/nrow(estimates_NLME)
  bias_perc_NLME[i]=bias_NLME[i]/true_value[i]
}
Coverage_NLME <- sapply(CI_NLME, mean)
result_NLME <- cbind(true_value,estimates = sapply(estimates_NLME, mean), se_m = sapply(SE_NLME, mean), se_s = sapply(estimates_NLME, sd), rMSE_NLME,bias_NLME,bias_perc_NLME, Coverage_NLME, NC = c(num_sim-nrow(estimates_NLME), NA, NA), time = c(sum(na.omit(time_NLME)), NA, NA))
rownames(result_NLME)=c("p1","lambda1","p2")

# Results table for LME4 
for (i in 1:num_sim){
  CI_LME4[i, ] = (estimates_LME4[i, ] - 1.96 * SE_LME4[i, ] <= true_value & estimates_LME4[i, ] + 1.96 * SE_LME4[i, ] >= true_value)
}
CI_LME4 <- na.omit(CI_LME4)
estimates_LME4 <- na.omit(estimates_LME4)
SE_LME4 <- na.omit(SE_LME4)
rMSE_LME4 <- bias_LME4 <- bias_perc_LME4 <- vector("double", length(true_value))
for (i in seq_along(true_value)){
  rMSE_LME4[i]=sqrt(sum((estimates_LME4[,i]-true_value[i])^2)/nrow(estimates_LME4))/true_value[i]
  bias_LME4[i]=sum(estimates_LME4[,i]-true_value[i])/nrow(estimates_LME4)
  bias_perc_LME4[i]=bias_LME4[i]/true_value[i]
}
Coverage_LME4 <- sapply(CI_LME4, mean)
result_LME4 <- cbind(true_value,estimates = sapply(estimates_LME4, mean), se_m = sapply(SE_LME4, mean), se_s = sapply(estimates_LME4, sd), rMSE_LME4,bias_LME4,bias_perc_LME4, Coverage_LME4, NC = c(num_sim-nrow(estimates_LME4), NA, NA), time = c(sum(na.omit(time_LME4)), NA, NA))
rownames(result_LME4)=c("p1","lambda1","p2")

# Latex Table
table_names <- c( "Method", "Time (s)", "NC","Parameter","True Value", "Estimate", "SE_M",  "SE_S", "Bias (%)","rMSE (%)","Coverage")
table <- data.frame(matrix(nrow=9, ncol=length(table_names)))
colnames(table) <- table_names
table$Parameter <- c("$P_1$","$\\lambda_1$","$P_2$","$P_1$","$\\lambda_1$","$P_2$","$P_1$","$\\lambda_1$","$P_2$")
table$`True Value` <- c(17, 4, 3, 17, 4, 3, 17, 4, 3)
table$Method <- c("SAEM",NA,NA,"NLME",NA,NA,"LME4",NA,NA)
table$Estimate <- c(result_SAEM[, 2], result_NLME[,2], result_LME4[,2])
table$SE_M <- c(result_SAEM[, 3], result_NLME[,3], result_LME4[,3])
table$SE_S <- c(result_SAEM[, 4], result_NLME[,4], result_LME4[,4])
table$`Bias (%)` <- c(result_SAEM[, 7]*100, result_NLME[,7]*100, result_LME4[,7]*100)
table$`rMSE (%)` <- c(result_SAEM[, 5]*100, result_NLME[,5]*100, result_LME4[,5]*100)
table$Coverage <-  c(result_SAEM[, 8]*100, result_NLME[,8]*100, result_LME4[,8]*100)
table$`Time (s)` <- c(result_SAEM[, 10], result_NLME[,10], result_LME4[,10])
table$NC <- c(result_SAEM[, 9], result_NLME[,9], result_LME4[,9])

latex_table<-xtable(table, type = "latex",align=c("cccccccccccc"))
digits(latex_table)<-c(0,2,1,0,1,1,2,2,2,2,2,2)
print(latex_table, file = "Table/Decay_n200_ni15.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,9))



