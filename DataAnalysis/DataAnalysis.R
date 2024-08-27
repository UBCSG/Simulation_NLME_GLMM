library(lixoftConnectors)
library(mvtnorm)
library(nlme)
library(MASS)
library(saemix)
library(tidyverse)
library(testit)
library(here)
library(docopt)
library(ggplot2)
library(lmtest)
library(gridExtra)
library(splines)
library(stats)
library(lme4)
library(GLDEX)
library(lme4)
library(xtable)
library(berryFunctions)
library(latex2exp)
initializeLixoftConnectors(software="monolix")
setwd(here::here("DataAnalysis"))

## Create results table 
decay_SAEM <- decay_NLME <- decay_LME4 <- data.frame(matrix(ncol = 5, nrow = 4))
colnames(decay_SAEM) <- colnames(decay_NLME) <- colnames(decay_LME4) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
decay_SAEM$Parameter <- decay_NLME$Parameter <- decay_LME4$Parameter <- c("$P_1$","$\\lambda_2$","\\beta", "$P_2$")

rebound_SAEM <- rebound_NLME <- rebound_LME4 <- data.frame(matrix(ncol = 5, nrow = 4))
colnames(rebound_SAEM) <- colnames(rebound_NLME) <- colnames(rebound_LME4) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
rebound_SAEM$Parameter <- rebound_NLME$Parameter <- rebound_LME4$Parameter <- c("$\\beta_1$","$\\beta_2$","$\\beta_3$","$\\beta_4$")

GLMM_SAEM <- GLMM_LME4 <- data.frame(matrix(ncol = 5, nrow = 4))
colnames(GLMM_SAEM) <-  colnames(GLMM_LME4) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
GLMM_SAEM$Parameter <- GLMM_LME4$Parameter <- c("$\\alpha_0$","$\\alpha_1$","$\\alpha_2$","\\alpha_3")

############ Decay model ==========
# SAEM 
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/decay2.csv'),
                      headerTypes =c("id","time","observation","contcov"),
                      observationTypes =list(log10 = "continuous"))
modelFile = paste0('Model/model_decay.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_new, modelFile = modelFile)
# set error model and observation distribution
setErrorModel(list(log10_ = "constant"))
setObservationDistribution(log10_ = "normal")
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setCorrelationBlocks(id = list(  c("p1", "b1", "p2") ) )
setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE)
setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal")
setCovariateModel(b1 = c(CD4 = TRUE))
setPopulationParameterInformation(p1_pop = list(initialValue = 18), 
                                  b1_pop = list(initialValue = 4), 
                                  p2_pop = list(initialValue = 3))
# run the estimation
setScenario(scenario)
runScenario()

report_SAEM = getSAEMiterations()$estimates
report_SAEM <- report_SAEM %>% 
  mutate(omega_p1 = omega_p1^2) %>% 
  mutate(omega_p2 = omega_p2^2) %>% 
  mutate(omega_b1 = omega_b1^2) 
report_SAEM <- rownames_to_column(report_SAEM, "IterationNum")
report_SAEM_long <- report_SAEM %>%
  pivot_longer(cols = 2:13, names_to = "Parameter", values_to = "Estimate")
report_SAEM_long$IterationNum <- as.numeric(report_SAEM_long$IterationNum)
report_SAEM_long$Parameter <- factor(report_SAEM_long$Parameter, levels=c("p1_pop","b1_pop", "beta_b1_CD4", "p2_pop", "omega_p1", "omega_p2","omega_b1", "corr_p1_b1","corr_p2_b1", "corr_p2_p1", "a", "convergenceIndicator"),
                                     labels=c(bquote(P[1]), bquote(lambda[1]), bquote(beta), bquote(P[2]), bquote(B[11]), bquote(B[22]), bquote(B[33]), bquote("cor(" ~ b[1][i]~","~b[3][i]~ ")"), bquote("cor("~b[2][i]~","~b[3][i]~")"), bquote("cor(" ~ b[1][i]~","~b[2][i]~ ")"), bquote(sigma),  bquote("CompleteLikelihood")))
ggplot(report_SAEM_long, aes(x = IterationNum, y = Estimate, group = Parameter)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", labeller = label_parsed) + # Create separate panels for each parameter
  labs(x = "Iteration", y = "Parameter Estimate", title = "Trajectories of Parameter Estimates") +
  theme_minimal()

# store the estimates
decay_SAEM$Estimates <- getEstimatedPopulationParameters()[1:4]
decay_SAEM$`Standard Error`<- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
decay_SAEM$`$z$-value` <- decay_SAEM$Estimates /decay_SAEM$`Standard Error`
decay_SAEM$`$p$-value` <- 2 * pnorm(abs(decay_SAEM$`$z$-value`[1:4]), lower.tail = FALSE)
sigma_SAEM <- getEstimatedPopulationParameters()[11]
B_SAEM <- diag(getEstimatedPopulationParameters()[5:7]^2)
B_SAEM[1,3] <- B_SAEM[3,1] <- getEstimatedPopulationParameters()[5] * getEstimatedPopulationParameters()[6] * getEstimatedPopulationParameters()[8]
B_SAEM[1,2] <- B_SAEM[2,1] <- getEstimatedPopulationParameters()[5] * getEstimatedPopulationParameters()[7] * getEstimatedPopulationParameters()[10]
B_SAEM[2,3] <- B_SAEM[3,2] <- getEstimatedPopulationParameters()[6] * getEstimatedPopulationParameters()[7] * getEstimatedPopulationParameters()[9]



# NLME
data.decay = read.csv("Data/Cleaned_Data/decay2.csv")
data.decay <- groupedData(log10 ~ time | PATIENT, data = data.decay)
nlme_control <- nlmeControl(msMaxIter = 50, msVerbose = TRUE)
nlme.decay <- nlme(log10 ~ log10(exp(p1-b1*time)+exp(p2)),
                   fixed = list(p1 ~ 1, b1 ~ CD4, p2 ~ 1),
                   random = p1+b1+p2 ~ 1,
                   data = data.decay,
                   start = c(17, 4, -0.1, 3),
                   control = nlme_control)


decay_NLME$Estimates  <- summary(nlme.decay)$tTable[,1]
decay_NLME$`Standard Error` <- summary(nlme.decay)$tTable[,2]
decay_NLME$`$z$-value` <- decay_NLME$Estimates/decay_NLME$`Standard Error`
decay_NLME$`$p$-value` <- 2 * pnorm(abs(decay_NLME$`$z$-value`), lower.tail = FALSE)
sigma_NLME <- VarCorr(nlme.decay)[4,2]
vcov_NLME <- as.data.frame(VarCorr(nlme.decay))
B_NLME <- diag(VarCorr(nlme.decay)[1:3, 1])
B_NLME[1, 3] <- B_NLME[3, 1] <- as.numeric(VarCorr(nlme.decay)[2, 3]) * as.numeric(VarCorr(nlme.decay)[1, 2]) * as.numeric(VarCorr(nlme.decay)[2, 2])
B_NLME[1, 2] <- B_NLME[2, 1] <- as.numeric(VarCorr(nlme.decay)[3, 3]) * as.numeric(VarCorr(nlme.decay)[1, 2]) * as.numeric(VarCorr(nlme.decay)[3, 2])
B_NLME[2, 3] <- B_NLME[3, 2] <- as.numeric(VarCorr(nlme.decay)[3, 4]) * as.numeric(VarCorr(nlme.decay)[3, 2]) * as.numeric(VarCorr(nlme.decay)[2, 2])
B_NLME

# u = evaluate({
#   nlme.decay <- nlme(log10 ~ log10(exp(p1-b1*time)+exp(p2)),
#                      fixed = list(p1 ~ 1, b1 ~ CD4, p2 ~ 1),
#                      random = p1+b1+p2 ~ 1,
#                      data = data.decay,
#                      start = c(17, 4, -0.1, 3),
#                      verbose = TRUE, control = nlmeControl(trace = TRUE))
# })
# 
# # get the output in character form
# output_str = u[[2]]
# 
# # check if everything worked
# cat(output_str)

# LME4
nform <- ~log10 (exp(p1-(b1+beta*CD4)*t)+ exp (p2))
nfun <- deriv(nform, namevec = c("p1", "b1","beta", "p2"), 
              function.arg = c("t","CD4", "p1", "b1", "beta", "p2"))
lme4.decay <- nlmer(log10 ~ nfun(time, CD4, p1, b1, beta, p2) ~ p1 + b1 + p2|PATIENT, 
                    data.decay, 
                    start = c(p1 = 18, b1 = 4, beta = -0.1, p2 = 3),
                    control = nlmerControl(optimizer = "nlminbwrap"))
decay_LME4$Estimates  <- summary(lme4.decay)$coefficients[, 1]
decay_LME4$`Standard Error` <- summary(lme4.decay)$coefficients[,2]
decay_LME4$`$z$-value` <- decay_LME4$Estimates/decay_LME4$`Standard Error`
decay_LME4$`$p$-value` <- 2 * pnorm(abs(decay_LME4$`$z$-value`), lower.tail = FALSE)
vcov_LME4 <- as.data.frame(VarCorr(lme4.decay))
sigma_LME4 <- vcov_LME4[7, 5]
B_LME4 <- diag(vcov_LME4[1:3, 4])
B_LME4[1, 2] <- B_LME4[2, 1] <- vcov_LME4[5, 4]
B_LME4[1, 3] <- B_LME4[3, 1] <- vcov_LME4[4, 4]
B_LME4[2, 3] <- B_LME4[3, 2] <- vcov_LME4[6, 4]
B_LME4

# Save results in LaTex
table <- rbind(decay_SAEM, decay_NLME, decay_LME4)
Method <- c("SAEM", NA, NA, NA, "NLME", NA, NA, NA, "LME4", NA, NA, NA) 
table <- cbind(Method, table)
latex_table<-xtable(table, type = "latex",align=c("ccccccc"))
digits(latex_table)<-c(0,0,2,2,2,2,3)
print(latex_table, file = "Results/Decay2.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,9))



# ############ Rebound model ==========
# # SAEM
# # Read in data for Monolix
# data_new = list(dataFile = paste0('Data/Cleaned_Data/rebound.csv'),
#                 headerTypes =c("id","time","observation"),
#                 observationTypes =list(log10 = "continuous"))
# modelFile = paste0('Model/model_rebound.txt')
# # create a new project by setting a data set and a structural model
# newProject(data = data_new, modelFile = modelFile)
# # set error model and observation distribution
# setErrorModel(list(log10_ = "constant"))
# setObservationDistribution(log10_= "normal")
# # set tasks in scenario
# scenario <- getScenario()
# scenario$tasks = c(populationParameterEstimation = T,
#                    conditionalModeEstimation = T,
#                    conditionalDistributionSampling = T,
#                    standardErrorEstimation=T,
#                    logLikelihoodEstimation=T)
# scenario$linearization = FALSE
# setIndividualParameterVariability(beta1 = TRUE,beta2 = TRUE,beta3 = TRUE,beta4 = TRUE)
# setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
# setPopulationParameterInformation(beta1_pop = list(initialValue = 3), 
#                                   beta2_pop = list(initialValue = 8), 
#                                   beta3_pop = list(initialValue = 2), 
#                                   beta4_pop = list(initialValue = 1.5))
# # run the estimation
# setScenario(scenario)
# runScenario()
# 
# # store the estimates
# rebound_SAEM$Estimates <- getEstimatedPopulationParameters()[1:4]
# rebound_SAEM$`Standard Error` <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
# rebound_SAEM$`$z$-value` <- rebound_SAEM$Estimates/rebound_SAEM$`Standard Error`
# rebound_SAEM$`$p$-value` <- 2 * pnorm(abs(rebound_SAEM$`$z$-value`), lower.tail = FALSE)
# 
# # NLME 
# data.rebound = read.csv("Data/Cleaned_Data/rebound.csv")
# data.rebound <- groupedData(log10 ~ time | PATIENT, data = data.rebound)
# nlme.rebound <- nlme(log10 ~ beta1 * time / (time  + exp( beta2 - beta3 * time)) + beta4,
#                      fixed = beta1+beta2+beta3+beta4~ 1,
#                      random = beta1+beta2+beta3+beta4 ~ 1,
#                      data = data.rebound,
#                      start = c(beta1 = 3, beta2 = 8, beta3 = 2, beta4 = 1.5))
# 
# rebound_NLME$Estimates  <- summary(nlme.rebound)$tTable[,1]
# rebound_NLME$`Standard Error` <- summary(nlme.rebound)$tTable[,2]
# rebound_NLME$`$z$-value` <- rebound_NLME$Estimates/rebound_NLME$`Standard Error`
# rebound_NLME$`$p$-value` <- 2 * pnorm(abs(rebound_NLME$`$z$-value`), lower.tail = FALSE)
# 
# # LME4
# nform <- ~ (beta1 * t / (t  + exp( beta2 - beta3 * t)) + beta4)
# nfun <- deriv(nform, namevec = c("beta1", "beta2", "beta3", "beta4"), 
#               function.arg = c("t", "beta1", "beta2", "beta3", "beta4"))
# lme4.rebound <- nlmer(log10 ~ nfun(time, beta1, beta2, beta3, beta4) ~ (beta1 + beta2 + beta3 + beta4| PATIENT), data.rebound, start = c(beta1 = 3, beta2 = 8, beta3 = 2, beta4 = 1.5))
# rebound_LME4$Estimates  <- summary(lme4.rebound)$coefficients[, 1]
# rebound_LME4$`Standard Error` <- summary(lme4.rebound)$coefficients[,2]
# rebound_LME4$`$z$-value` <- rebound_LME4$Estimates/rebound_LME4$`Standard Error`
# rebound_LME4$`$p$-value` <- 2 * pnorm(abs(rebound_LME4$`$z$-value`), lower.tail = FALSE)
# 
# # Save results in LaTex
# table <- rbind(rebound_SAEM, rebound_NLME, rebound_LME4)
# Method <- c("SAEM", NA, NA, NA, "NLME", NA, NA, NA, "LME4", NA, NA, NA) 
# table <- cbind(Method, table)
# latex_table<-xtable(table, type = "latex",align=c("ccccccc"))
# digits(latex_table)<-c(0,0,2,2,2,2,3)
# print(latex_table, file = "Results/rebound.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,9))

############ GLMM ==========
# Choose the starting values
data.GLMM = read.csv("Data/Cleaned_Data/GLMM2.csv")
data.GLMM <- groupedData(censor ~ time | PATIENT, data = data.GLMM)
glm <- glm(censor ~ time + CD4_st + I(time^2),
  data = data.GLMM, family = binomial)
summary(glm)
# SAEM
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/GLMM2.csv'),
                headerTypes =c("id","observation","regressor","time"),
                observationTypes = list(censor = "categorical"))
modelFile = paste0('Model/model_GLMM2.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_new, modelFile = modelFile)
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setCorrelationBlocks(id = list(  c("gamma0", "gamma1", "gamma2", "gamma3") ) )
setIndividualParameterVariability(gamma0 = TRUE, gamma1 = TRUE, gamma2 = TRUE, gamma3 = TRUE)
setIndividualParameterDistribution(gamma0="normal", gamma1="normal", gamma2="normal", gamma3 = "normal")
setPopulationParameterInformation(gamma0_pop = list(initialValue = -1), 
                                  gamma1_pop  = list(initialValue = 0.2), 
                                  gamma2_pop  = list(initialValue = 0.4),
                                  gamma3_pop = list(initialValue = -0.01))
# run the estimation
setScenario(scenario)
runScenario()

# store the estimates
GLMM_SAEM$Estimates <- getEstimatedPopulationParameters()[1:4]
GLMM_SAEM$`Standard Error` <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
GLMM_SAEM$`$z$-value` <- GLMM_SAEM$Estimates/GLMM_SAEM$`Standard Error`
GLMM_SAEM$`$p$-value` <- 2 * pnorm(abs(GLMM_SAEM$`$z$-value`), lower.tail = FALSE)

A_SAEM <- diag(getEstimatedPopulationParameters()[5:8]^2)
A_SAEM[1,3] <- A_SAEM[3,1] <- getEstimatedPopulationParameters()[5] * getEstimatedPopulationParameters()[7] * getEstimatedPopulationParameters()[10]
A_SAEM[1,2] <- A_SAEM[2,1] <- getEstimatedPopulationParameters()[5] * getEstimatedPopulationParameters()[6] * getEstimatedPopulationParameters()[9]
A_SAEM[1,4] <- A_SAEM[4,1] <- getEstimatedPopulationParameters()[5] * getEstimatedPopulationParameters()[8] * getEstimatedPopulationParameters()[11]
A_SAEM[2,3] <- A_SAEM[3,2] <- getEstimatedPopulationParameters()[6] * getEstimatedPopulationParameters()[7] * getEstimatedPopulationParameters()[12]
A_SAEM[2,4] <- A_SAEM[4,2] <- getEstimatedPopulationParameters()[6] * getEstimatedPopulationParameters()[8] * getEstimatedPopulationParameters()[13]
A_SAEM[3,4] <- A_SAEM[4,3] <- getEstimatedPopulationParameters()[7] * getEstimatedPopulationParameters()[8] * getEstimatedPopulationParameters()[14]

report_SAEM = getSAEMiterations()$estimates
report_SAEM <- report_SAEM %>% 
  mutate(omega_gamma0 = omega_gamma0^2) %>% 
  mutate(omega_gamma1 = omega_gamma1^2) %>% 
  mutate(omega_gamma2 = omega_gamma2^2) %>% 
  mutate(omega_gamma3 = omega_gamma3^2)
report_SAEM <- rownames_to_column(report_SAEM, "IterationNum")
report_SAEM_long <- report_SAEM %>%
  pivot_longer(cols = 2:16, names_to = "Parameter", values_to = "Estimate")
report_SAEM_long$IterationNum <- as.numeric(report_SAEM_long$IterationNum)
report_SAEM_long$Parameter <- factor(report_SAEM_long$Parameter, levels=c("gamma0_pop","gamma1_pop", "gamma2_pop", "gamma3_pop", "omega_gamma0", "omega_gamma1","omega_gamma2", "omega_gamma3","corr_gamma1_gamma0", "corr_gamma2_gamma0", "corr_gamma3_gamma0", "corr_gamma2_gamma1", "corr_gamma3_gamma1", "corr_gamma3_gamma2", "convergenceIndicator"),
                                     labels=c(bquote(alpha[0]),bquote(alpha[1]),bquote(alpha[2]), bquote(alpha[3]), bquote(A[11]), bquote(A[22]), bquote(A[33]),bquote(A[44]), bquote("cor(" ~ a[0][i]~","~a[1][i]~ ")"), bquote("cor("~a[0][i]~","~a[2][i]~")"), bquote("cor(" ~ a[0][i]~","~a[3][i]~ ")"), bquote("cor(" ~ a[1][i]~","~a[2][i]~ ")"), bquote("cor(" ~ a[1][i]~","~a[3][i]~ ")"), bquote("cor(" ~ a[2][i]~","~a[3][i]~ ")"),   bquote("CompleteLikelihood")))
ggplot(report_SAEM_long, aes(x = IterationNum, y = Estimate, group = Parameter)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", labeller = label_parsed) + # Create separate panels for each parameter
  labs(x = "Iteration", y = "Parameter Estimate", title = "Trajectories of Parameter Estimates") +
  theme_minimal()

# # SAEM
# # Read in data for Monolix
# data_new = list(dataFile = paste0('Data/Cleaned_Data/GLMM.csv'),
#                 headerTypes =c("id","observation","regressor","time"),
#                 observationTypes = list(CD4 = "categorical"))
# modelFile = paste0('Model/model_GLMM3.txt')
# # create a new project by setting a data set and a structural model
# newProject(data = data_new, modelFile = modelFile)
# # set tasks in scenario
# scenario <- getScenario()
# scenario$tasks = c(populationParameterEstimation = T,
#                    conditionalModeEstimation = T,
#                    conditionalDistributionSampling = T,
#                    standardErrorEstimation=T,
#                    logLikelihoodEstimation=T)
# scenario$linearization = FALSE
# setIndividualParameterVariability(gamma0 = TRUE, gamma1 = TRUE, gamma2 = TRUE, gamma3 = TRUE)
# setIndividualParameterDistribution(gamma0="normal", gamma1="normal", gamma2="normal", gamma3 = "normal")
# setPopulationParameterInformation(gamma0_pop = list(initialValue = -10), 
#                                   gamma1_pop  = list(initialValue = 3), 
#                                   gamma2_pop  = list(initialValue = 0.1),
#                                   gamma3_pop = list(initialValue = -0.1))
# # setCovariateModel(gamma0 = c(viral_st = TRUE))
# # run the estimation
# setScenario(scenario)
# runScenario()
# 
# # store the estimates
# GLMM_SAEM$Estimates <- getEstimatedPopulationParameters()[1:4]
# GLMM_SAEM$`Standard Error` <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
# GLMM_SAEM$`$z$-value` <- GLMM_SAEM$Estimates/GLMM_SAEM$`Standard Error`
# GLMM_SAEM$`$p$-value` <- 2 * pnorm(abs(GLMM_SAEM$`$z$-value`), lower.tail = FALSE)

# LME4
lme4.GLMM <- glmer(
  censor ~ time + CD4_st + I(time^2) + ( 1 + time + CD4_st  + I(time^2) |PATIENT),
  data = data.GLMM, family = binomial, start = list(fixef = c(-1, 0.2, 0.4, -0.01)))
GLMM_LME4$Estimates  <- summary(lme4.GLMM)$coefficients[, 1]
GLMM_LME4$`Standard Error` <- summary(lme4.GLMM)$coefficients[,2]
GLMM_LME4$`$z$-value` <- GLMM_LME4$Estimates/GLMM_LME4$`Standard Error`
GLMM_LME4$`$p$-value` <- 2 * pnorm(abs(GLMM_LME4$`$z$-value`), lower.tail = FALSE)
vcov_LME4 <- as.data.frame(VarCorr(lme4.GLMM))
A_LME4 <- diag(vcov_LME4[1:4, 4])
A_LME4[1, 2] <- A_LME4[2, 1] <- vcov_LME4[5, 4]
A_LME4[1, 3] <- A_LME4[3, 1] <- vcov_LME4[6, 4]
A_LME4[1, 4] <- A_LME4[4, 1] <- vcov_LME4[7, 4]
A_LME4[2, 3] <- A_LME4[3, 2] <- vcov_LME4[8, 4]
A_LME4[2, 4] <- A_LME4[4, 2] <- vcov_LME4[9, 4]
A_LME4[4, 3] <- A_LME4[3, 4] <- vcov_LME4[10, 4]
A_LME4

# data.GLMM = read.csv("Data/Cleaned_Data/GLMM.csv")
# data.GLMM <- groupedData(CD4 ~ time | PATIENT, data = data.GLMM)
# lme4.GLMM <- glmer(
#   CD4 ~ time + viral_st + I(time^2) + (1 +  time + viral_st + I(time^2)|PATIENT),
#   data = data.GLMM, family = binomial)
# GLMM_LME4$Estimates  <- summary(lme4.GLMM)$coefficients[, 1]
# GLMM_LME4$`Standard Error` <- summary(lme4.GLMM)$coefficients[,2]
# GLMM_LME4$`$z$-value` <- GLMM_LME4$Estimates/GLMM_LME4$`Standard Error`
# GLMM_LME4$`$p$-value` <- 2 * pnorm(abs(GLMM_LME4$`$z$-value`), lower.tail = FALSE)

# Save results in LaTex
table <- rbind(GLMM_SAEM, GLMM_LME4)
Method <- c("SAEM", NA, NA,NA, "LME4", NA, NA,NA) 
table <- cbind(Method, table)
latex_table<-xtable(table, type = "latex",align=c("ccccccc"))
digits(latex_table)<-c(0,0,3,3,3,3,3)
print(latex_table, file = "Results/GLMM2.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,6))
