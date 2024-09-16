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
decay_SAEM <- decay_NLME <- decay_LME4 <- data.frame(matrix(ncol = 5, nrow = 8))
colnames(decay_SAEM) <- colnames(decay_NLME) <- colnames(decay_LME4) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
decay_SAEM$Parameter <- decay_NLME$Parameter <- decay_LME4$Parameter <- c("$P_1$","$\\lambda_2$","$\\beta$", "$P_2$", "$\\sigma$", "$B_{11}$", "$B_{22}$", "$B_{33}$")

GLMM_SAEM <- GLMM_LME4 <- data.frame(matrix(ncol = 5, nrow = 8))
colnames(GLMM_SAEM) <-  colnames(GLMM_LME4) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
GLMM_SAEM$Parameter <- GLMM_LME4$Parameter <- c("$\\alpha_0$","$\\alpha_1$","$\\alpha_2$","$\\alpha_3$", "$A_{11}$", "$A_{22}$", "$A_{33}$", "$A_{44}$")

############ NLME: Viral decay model ==========
# y_ij = log10(exp(P_{1i}-lambda_{1i}*t_{ij}) + exp(P_{2i})) + e_{ij}
# P_{1i} = P_1 + b_{1i}; P_{2i} = P_2 + b_{2i}; lambda_{1i} = lambda_1 + beta * CD4_i + b_{3i},
# b_i ~ N(0, B), e_{ij} ~ N(0, sigma^2)

# SAEM 
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/decay.csv'),
                      headerTypes =c("id","time","observation","contcov"),
                      observationTypes =list(log10 = "continuous"))
modelFile = paste0('Model/model_decay.txt')
# Create a new project by setting a data set and a structural model
newProject(data = data_new, modelFile = modelFile)
# Set error model and observation distribution
setErrorModel(list(log10_ = "constant"))
setObservationDistribution(log10_ = "normal")
# Set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setCorrelationBlocks(id = list(c("p1", "b1", "p2")))
setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE)
setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal")
setCovariateModel(b1 = c(CD4 = TRUE))
setPopulationParameterInformation(p1_pop = list(initialValue = 18), 
                                  b1_pop = list(initialValue = 4), 
                                  p2_pop = list(initialValue = 3))
# Run the estimation
setScenario(scenario)
runScenario()

# Draw convergence plot
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

# Store the estimates
decay_SAEM$Estimates[1:4] <- getEstimatedPopulationParameters()[1:4]
decay_SAEM$`Standard Error`[1:4] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
decay_SAEM$`$z$-value`[1:4] <- decay_SAEM$Estimates[1:4] /decay_SAEM$`Standard Error`[1:4]
decay_SAEM$`$p$-value`[1:4] <- 2 * pnorm(abs(decay_SAEM$`$z$-value`[1:4]), lower.tail = FALSE)
decay_SAEM$Estimates[5] <- getEstimatedPopulationParameters()[11]
decay_SAEM$Estimates[6:8] <- getEstimatedPopulationParameters()[c(5, 7, 6)]^2

# NLME
data.decay = read.csv("Data/Cleaned_Data/decay.csv")
data.decay <- groupedData(log10 ~ time | PATIENT, data = data.decay)
nlme_control <- nlmeControl(msMaxIter = 50, msVerbose = TRUE)
nlme.decay <- nlme(log10 ~ log10(exp(p1-b1*time)+exp(p2)),
                   fixed = list(p1 ~ 1, b1 ~ CD4, p2 ~ 1),
                   random = p1+b1+p2 ~ 1,
                   data = data.decay,
                   start = c(17, 4, -0.1, 3),
                   control = nlme_control)

# Store the estimates
decay_NLME$Estimates[1:4]  <- summary(nlme.decay)$tTable[,1]
decay_NLME$`Standard Error`[1:4] <- summary(nlme.decay)$tTable[,2]
decay_NLME$`$z$-value`[1:4] <- decay_NLME$Estimates[1:4]/decay_NLME$`Standard Error`[1:4]
decay_NLME$`$p$-value`[1:4] <- 2 * pnorm(abs(decay_NLME$`$z$-value`[1:4]), lower.tail = FALSE)
decay_NLME$Estimates[5] <- VarCorr(nlme.decay)[4, 2]
decay_NLME$Estimates[6:8] <- VarCorr(nlme.decay)[c(1, 3, 2), 1]

# LME4
nform <- ~log10 (exp(p1-(b1+beta*CD4)*t)+ exp (p2))
nfun <- deriv(nform, namevec = c("p1", "b1","beta", "p2"), 
              function.arg = c("t","CD4", "p1", "b1", "beta", "p2"))
lme4.decay <- nlmer(log10 ~ nfun(time, CD4, p1, b1, beta, p2) ~ p1 + b1 + p2|PATIENT, 
                    data.decay, 
                    start = c(p1 = 18, b1 = 4, beta = -0.1, p2 = 3),
                    control = nlmerControl(optimizer = "nlminbwrap"))

# Store the estimates
decay_LME4$Estimates[1:4]  <- summary(lme4.decay)$coefficients[, 1]
decay_LME4$`Standard Error`[1:4] <- summary(lme4.decay)$coefficients[,2]
decay_LME4$`$z$-value`[1:4] <- decay_LME4$Estimates[1:4]/decay_LME4$`Standard Error`[1:4]
decay_LME4$`$p$-value`[1:4] <- 2 * pnorm(abs(decay_LME4$`$z$-value`[1:4]), lower.tail = FALSE)
vcov_LME4 <- as.data.frame(VarCorr(lme4.decay))
decay_LME4$Estimates[5] <- vcov_LME4[7, 5]
decay_LME4$Estimates[6:8] <- vcov_LME4[c(1, 3, 2), 4]


# Save results in LaTex
table <- rbind(decay_SAEM, decay_NLME, decay_LME4)
Method <- c("SAEM", rep(NA, 7), "NLME", rep(NA, 7), "LME4", rep(NA, 7)) 
table$Estimates <- as.numeric(table$Estimates)
table <- cbind(Method, table)
latex_table <- xtable(table, type = "latex", align=c("ccccccc"))
digits(latex_table) <- c(0,2,2,2,2,2,3)
print(latex_table, file = "Results/Decay.tex",include.rownames=FALSE, sanitize.text.function = function(x){x}, hline.after = c(-1,0,24))



############ GLMM ==========
# logit(P(z_{ij} = 1)) = (alpha_0 + a_{0i}) + (alpha_1 + a_{1i}) * t_{ij} + (alpha_2 + a_{2i}) * w_{ij} + (alpha_3 + a_{3i}) * t_{ij}^2
# z_{ij} = 0 if the viral load is observed and z_{ij} = 1 if the viral load is censored.
# w_{ij} is the standardied CD4 value at time t_{ij} and a_i ~ N(0, A).

# Choose the starting values
data.GLMM = read.csv("Data/Cleaned_Data/GLMM.csv")
data.GLMM <- groupedData(censor ~ time_raw | PATIENT, data = data.GLMM)
glm <- glm(censor ~ time_raw + CD4_st + I(time_raw^2),
  data = data.GLMM, family = binomial)
summary(glm)
# SAEM
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/GLMM.csv'),
                headerTypes =c("id","time","observation","regressor"),
                observationTypes = list(censor = "categorical"))
modelFile = paste0('Model/model_GLMM.txt')
# Create a new project by setting a data set and a structural model
newProject(data = data_new, modelFile = modelFile)
# Set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setCorrelationBlocks(id = list( c("alpha0", "alpha1", "alpha2", "alpha3") ) )
setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, alpha3 = TRUE)
setIndividualParameterDistribution(alpha0="normal", alpha1="normal", alpha2="normal", alpha3 = "normal")
setPopulationParameterInformation(alpha0_pop = list(initialValue = -1), 
                                  alpha1_pop  = list(initialValue = 0.2), 
                                  alpha2_pop  = list(initialValue = 0.4),
                                  alpha3_pop = list(initialValue = -0.01))
# Run the estimation
setScenario(scenario)
runScenario()

# Store the estimates
GLMM_SAEM$Estimates[1:4] <- getEstimatedPopulationParameters()[1:4]
GLMM_SAEM$`Standard Error`[1:4] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:4])
GLMM_SAEM$`$z$-value`[1:4] <- GLMM_SAEM$Estimates[1:4]/GLMM_SAEM$`Standard Error`[1:4]
GLMM_SAEM$`$p$-value`[1:4] <- 2 * pnorm(abs(GLMM_SAEM$`$z$-value`[1:4]), lower.tail = FALSE)
GLMM_SAEM$Estimates[5:8] <- getEstimatedPopulationParameters()[5:8]^2

# Draw the convergence plot
report_SAEM = getSAEMiterations()$estimates
report_SAEM <- report_SAEM %>%
  mutate(omega_alpha0 = omega_alpha0 ^ 2) %>%
  mutate(omega_alpha1 = omega_alpha1 ^ 2) %>%
  mutate(omega_alpha2 = omega_alpha2 ^ 2) %>%
  mutate(omega_alpha3 = omega_alpha3 ^ 2)
report_SAEM <- rownames_to_column(report_SAEM, "IterationNum")
report_SAEM_long <- report_SAEM %>%
  pivot_longer(cols = 2:16,
               names_to = "Parameter",
               values_to = "Estimate")
report_SAEM_long$IterationNum <- as.numeric(report_SAEM_long$IterationNum)
report_SAEM_long$Parameter <- factor(report_SAEM_long$Parameter,
  levels = c(
    "alpha0_pop",
    "alpha1_pop",
    "alpha2_pop",
    "alpha3_pop",
    "omega_alpha0",
    "omega_alpha1",
    "omega_alpha2",
    "omega_alpha3",
    "corr_alpha1_alpha0",
    "corr_alpha2_alpha0",
    "corr_alpha3_alpha0",
    "corr_alpha2_alpha1",
    "corr_alpha3_alpha1",
    "corr_alpha3_alpha2",
    "convergenceIndicator"
  ),
  labels = c(
    bquote(alpha[0]),
    bquote(alpha[1]),
    bquote(alpha[2]),
    bquote(alpha[3]),
    bquote(A[11]),
    bquote(A[22]),
    bquote(A[33]),
    bquote(A[44]),
    bquote("cor(" ~ a[0][i] ~ "," ~ a[1][i] ~ ")"),
    bquote("cor(" ~ a[0][i] ~ "," ~ a[2][i] ~ ")"),
    bquote("cor(" ~ a[0][i] ~ "," ~ a[3][i] ~ ")"),
    bquote("cor(" ~ a[1][i] ~ "," ~ a[2][i] ~ ")"),
    bquote("cor(" ~ a[1][i] ~ "," ~ a[3][i] ~ ")"),
    bquote("cor(" ~ a[2][i] ~ "," ~ a[3][i] ~ ")"),
    bquote("CompleteLikelihood")
  )
)
ggplot(report_SAEM_long, aes(x = IterationNum, y = Estimate, group = Parameter)) +
  geom_line() +
  facet_wrap( ~ Parameter, scales = "free_y", labeller = label_parsed) + # Create separate panels for each parameter
  labs(x = "Iteration", y = "Parameter Estimate", title = "Trajectories of Parameter Estimates") +
  theme_minimal()


# LME4
lme4.GLMM <- glmer(
  censor ~ time_raw + CD4_st + I(time_raw^2) + ( 1 + time_raw + CD4_st  + I(time_raw^2) |PATIENT),
  data = data.GLMM, family = binomial, start = list(fixef = c(-1, 0.2, 0.4, -0.01)))
GLMM_LME4$Estimates[1:4]  <- summary(lme4.GLMM)$coefficients[, 1]
GLMM_LME4$`Standard Error`[1:4] <- summary(lme4.GLMM)$coefficients[,2]
GLMM_LME4$`$z$-value`[1:4] <- GLMM_LME4$Estimates[1:4]/GLMM_LME4$`Standard Error`[1:4]
GLMM_LME4$`$p$-value`[1:4] <- 2 * pnorm(abs(GLMM_LME4$`$z$-value`[1:4]), lower.tail = FALSE)
vcov_LME4 <- as.data.frame(VarCorr(lme4.GLMM))
GLMM_LME4$Estimates[5:8] <- vcov_LME4[1:4, 4]



# Save results in LaTex
table <- rbind(GLMM_SAEM, GLMM_LME4)
Method <- c("SAEM", rep(NA, 7), "LME4", rep(NA, 7)) 
table <- cbind(Method, table)
latex_table<-xtable(table, type = "latex",align=c("ccccccc"))
digits(latex_table)<-c(0,0,3,3,3,3,3)
print(latex_table, file = "Results/GLMM.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,16))
