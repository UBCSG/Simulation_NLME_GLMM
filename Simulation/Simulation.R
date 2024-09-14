rm(list=ls())

#======================== source all functions =====================
all.files = list.files(path = here::here("Simulation/src"), pattern="*.R$")
file.sources <- paste0(here::here("Simulation/src"), "/", all.files)
sapply(file.sources, source)

#======================== Simulation of Decay Model ====================
# Decay model: y_ij = log10(exp(p_1i - lambda_1i * t_ij) + exp(p_2i)) + e_ij
# Random effects: p_1i = p1 + b_1i, p_2i = p1 + b_2i, lambda_1i = lambda1 + b_3i,
# where e_ij ~ N(0, sigma1^2) and b_i ~ N(0, D) with diagonal elements being D11, D22, and D33.
# Default values: p1 = 17, lambda1 = 4, p2 = 3, D = diag(c(2, 0.1, 0.1)), sigma1 = 0.1, 
# time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 20.
SimFunDecay(p1 = 17, lambda1 = 4, p2 = 3, D = diag(c(2, 0.1, 0.1)),
            sigma1 = 0.1, time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 20)


##=====================Simulation of Rebound Model=======================
# Rebound model: y_ij = beta_1i * t_ij / (t_ij  + exp( beta_2i - beta_3i * t_ij)) + beta_4i / (1 + exp( beta_5i * t)) + e_ij
# Random effects: beta_ki = betak + tau_ki, k = 1, 2, 3, 4, 5
# where e_ij ~ N(0, sigma2^2) and tau_i ~ N(0, B) with diagonal elements being B11, B22, B33, and B44.
# Default values: beta1 = 4, beta2 = 5.7, beta3 = 2.1, beta4 = 1.9, beta5 = 0.4, 
# B = diag(c(0.2, 0, 0.1, 0, 0)), sigma = 0.1, time = c(0.5, 2.9, 4.8, 7, 10.2), num_sim = 500, num_patient = 50.
SimFunRebound(beta1 = 4, beta2 = 5.7, beta3 = 2.1, beta4 = 1.9, beta5 = 0.4, 
               B = diag(c(0.2, 0.1, 0.1, 0, 0)), sigma = 0.1, 
               time = c(0.5, 2.9, 4.8, 7, 10.2), num_sim = 50, num_patient = 50)



