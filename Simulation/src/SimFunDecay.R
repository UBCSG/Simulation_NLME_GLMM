##=====================Decay Model=======================
# Decay model: y_ij = log10(exp(p_1i - lambda_1i * t_ij) + exp(p_2i)) + e_ij
# Random effects: p_1i = p1 + b_1i, p_2i = p1 + b_2i, lambda_1i = lambda1 + b_3i,
# where e_ij ~ N(0, sigma1^2) and b_i ~ N(0, D) with diagonal elements being D11, D22, and D33.

##=====================Setting======================
SimFunDecay <- function(p1 = 17, lambda1 = 4, p2 = 3, D = diag(c(2, 0.1, 0.1)),
                        sigma = 0.1, time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 20) {
  # Find which parameter has a random effect
  D_diag <- which(diag(D) != 0)
  D_name <- ""
  for (i in 1:length(D_diag)){
    D_name = paste(D_name, "_D", D_diag[i], D_diag[i], "_", diag(D)[i], sep = "")
  }
  
  # Create a folder to store simulation results
  folder.name = folder.name = paste("Decay_n", num_patient, "_ni", length(time), "_p1_", p1, "_lambda1_", lambda1,
                                    "_p2_", p2, "_sigma", sigma, D_name, "_numsim", num_sim, sep="")
  dir.create(file.path("Simulation/", folder.name), showWarnings = FALSE)
  
  # Parameter estimates
  estimates_SAEM <- estimates_NLME <- estimates_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(estimates_SAEM) <- colnames(estimates_NLME) <- colnames(estimates_LME4) <- c("p1", "lambda1", "p2")
  # Variance of random effects
  D_SAEM <- D_NLME <- D_LME4 <- data.frame(matrix(nrow = num_sim, ncol = length(D_diag)))
  colnames(D_SAEM) <- colnames(D_NLME) <- colnames(D_LME4) <- c("D11", "D22", "D33")
  # SE of the parameter estimates
  SE_SAEM <- SE_NLME <- SE_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(SE_SAEM) <- colnames(SE_NLME) <- colnames(SE_LME4) <- c("p1", "lambda1", "p2")
  # Residual SD (sigma_1)
  sigma_SAEM <- sigma_NLME <- sigma_LME4 <- c(rep(NA, num_sim))
  # Computation time
  time_SAEM <- time_NLME <- time_LME4 <- c(rep(NA, num_sim))

  for (i in 1:num_sim) {
    # Print the iteration number
    print(glue::glue("i=", i))
    
    # Simulate dataset for viral load before ART interruption
    # Create a vector 'PATIENT' where each patient ID is repeated for each time point
    PATIENT <- rep(1:num_patient, each = length(time))
    # Create a vector 'day' that repeats the time points for each patient
    day <- rep(time, num_patient)
    # Combine the 'PATIENT' and 'day' vectors into a matrix 'data'
    data <- cbind(PATIENT, day)
    
    # Simulate b_i and e_ij
    bi_sim <- mvrnorm(num_patient, c(0, 0, 0), D) # bi ~ N(0, D)
    colnames(bi_sim) = c("b1_sim", "b2_sim", "b3_sim") # p_1i = p1 + b_1i, p_2i = p1 + b_2i, lambda_1i = lambda1 + b_3i
    bi_sim <- cbind(PATIENT = c(1:num_patient), bi_sim)
    data <- merge(data, bi_sim, by = "PATIENT")
    data <- data %>%
      mutate(e = rnorm(nrow(data), 0, sigma)) # e_ij ~ N(0, sigma^2)
    
    # Simulate viral load based on y_ij = log10(exp(p_1i - lambda_1i * t_ij) + exp(p_2i)) + e_ij
    data_sim <- data %>%
      mutate(viral = log10(exp(p1 + b1_sim - (lambda1 + b3_sim) * day) + exp(p2 + b2_sim)) + e) %>%
      dplyr::select(-b1_sim, -b2_sim, -b3_sim, -e)

    # Plot of simulated viral load 
    ggplot(data_sim[PATIENT %in% c(1:5), ], aes(x = day, y = viral)) +
      geom_point() +
      geom_line(aes(group = PATIENT)) +
      scale_x_continuous("Time") +
      scale_y_continuous(bquote("Viral load (in" ~ log[10] ~ "-scale)")) +
      theme_classic() + theme(
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90")
      )
    
    # Save the simulated data locally
    write.table(data_sim, paste0("Simulation/", folder.name, "/data.txt"), sep = "," , quote = FALSE, row.names = FALSE)
    
    # ================ Fit viral decay model using SAEM =================
    start.time <- Sys.time()
    data = list(
      dataFile = paste0("Simulation/", folder.name, "/data.txt"),
      headerTypes = c("id", "time", "observation"),
      observationTypes = list(viral = "continuous")
    )
    modelFile = paste0('Model/model_decay.txt')
    newProject(data = data, modelFile = modelFile)
    # Set the error model type to be used for the observation model
    setErrorModel(viral = "constant") # obs = pred + a*err, error ~ N(0, 1)
    # Set observation model distribution
    setObservationDistribution(viral = "normal") 
    # Set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(
      populationParameterEstimation = T,
      conditionalModeEstimation = T,
      conditionalDistributionSampling = T,
      standardErrorEstimation = T,
      logLikelihoodEstimation = T
    )
    scenario$linearization = FALSE
    setIndividualParameterVariability(p1 = TRUE, lambda1 = TRUE, p2 = TRUE)
    setPopulationParameterInformation(
      p1_pop = list(initialValue = p1),
      lambda1_pop = list(initialValue = lambda1),
      p2_pop = list(initialValue = p2)
    )
    setIndividualParameterDistribution(lambda1 = "normal", p1 = "normal", p2 = "normal")
    setScenario(scenario)
    # Run the estimation
    runScenario()
    
    # Store the estimates in table
    estimates_SAEM[i, ] <- getEstimatedPopulationParameters()[1:3]
    SE_SAEM[i, ] <- getEstimatedStandardErrors()$stochasticApproximation[1:3, 2]
    D_SAEM[i, ] <- getEstimatedPopulationParameters()[c(4, 6, 5)] ^ 2
    sigma_SAEM[i] <- getEstimatedPopulationParameters()[7]
    end.time <- Sys.time()
    time_SAEM[i] <- end.time - start.time
    
    # ================ Fit viral decay model using NLME =================
    start.time <- Sys.time()
    data_nlme <- groupedData(viral ~ day | PATIENT, data = data_sim)
  
    tryCatch({
      model_nlme <- nlme(
        viral ~ log10(exp(p1 - lambda1 * day) + exp(p2)),
        fixed = p1 + lambda1 + p2 ~ 1,
        random = pdDiag(list(p1 ~ 1, lambda1 ~ 1, p2 ~ 1)),
        data = data_nlme,
        start = c(p1 = p1, b1 = lambda1, p2 = p2)
      )
      
      estimates_NLME[i, ] <- summary(model_nlme)$tTable[, 1]
      SE_NLME[i, ] <- summary(model_nlme)$tTable[, 2]
      D_NLME[i, ] <- as.numeric(VarCorr(model_nlme)[1:3, 1])
      sigma_NLME[i] <- as.numeric(VarCorr(model_nlme)[4, 2])
    }, error = function(e) {
      estimates_NLME[i, ] <<- SE_NLME[i, ] <<- D_NLME[i, ] <<- sigma_NLME[i]  <<- NA
    })
    
    end.time <- Sys.time()
    time_NLME[i] <- end.time - start.time
    
    # ================ Fit viral decay model using LME4 =================
    start.time <- Sys.time()
    nform <- ~ log10 (exp(p1 - lambda1 * t) + exp (p2))
    nfun <- deriv(
      nform,
      namevec = c("p1", "lambda1", "p2"),
      function.arg = c("t", "p1", "lambda1", "p2")
    )
    tryCatch({
    model_lme4 <- nlmer(
      viral ~ nfun(day, p1, lambda1, p2) ~ (p1 | PATIENT) + (p2 | PATIENT) + (lambda1 | PATIENT),
      data_sim,
      start = c(p1 = p1, lambda1 = lambda1, p2 = p2)
    )
    estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
    SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
    D_LME4[i, ] <- as.numeric(as.data.frame(VarCorr(model_lme4))[1:3, 4])
    sigma_LME4[i] <- as.numeric(as.data.frame(VarCorr(model_lme4))[4, 5])
    }, error = function(e) {estimates_LME4[i, ] <<- SE_LME4[i, ] <<- D_LME4[i, ] <<- sigma_LME4[i] <<- NA})
    end.time <- Sys.time()
    time_LME4[i] <- end.time - start.time
  }
  
  # Create a folder to store the simulation results
  dir.create(file.path("Simulation/", folder.name, "Results"), showWarnings = FALSE)
  SE_SAEM[SE_SAEM == "NaN"] <- NA
  SE_SAEM <- mutate_all(SE_SAEM, function(x) as.numeric(as.character(x)))
  files_est <- list(estimates_SAEM = estimates_SAEM, estimates_NLME = estimates_NLME, estimates_LME4 = estimates_LME4,
                SE_SAEM = SE_SAEM, SE_NLME = SE_NLME, SE_LME4 = SE_LME4,
                D_SAEM = D_SAEM, D_NLME = D_NLME, D_LME4 = D_LME4)
  files_sigma <- list(sigma_SAEM = sigma_SAEM, sigma_NLME = sigma_NLME, sigma_LME4 = sigma_LME4)
  files_time <- list(time_SAEM = time_SAEM, time_NLME = time_NLME, time_LME4 = time_LME4)
  # Remove NA's
  files_est <- lapply(files_est, function(df) df[!apply(is.na(df), 1, all), ])
  files_sigma <- lapply(files_sigma, na.omit)
  files <- c(files_est, files_sigma, files_time)
  # Save simulation results
  for (i in seq_along(files)) {
    write.table(files[[i]], file = paste0("Simulation/", folder.name, "/Results/", names(files)[i], ".csv"))
  }
  
  #======================Calculate MSE and bias===================
  true_value <- c(p1, lambda1, p2)
  all_parameters <- c(p1, lambda1, p2, D[1, 1], D[2, 2], D[3, 3], sigma)

  # Results table for SAEM
  estimates_SAEM_all <- cbind(estimates_SAEM, D_SAEM, sigma_SAEM)
  # In case of any non-convergence, number of valid simulation iteration is
  num_SAEM <- nrow(estimates_SAEM)
  # CI_SAEM is a data frame containing boolean, TRUE is the estimate is within 95% CI, FALSE otherwise.
  CI_SAEM <- data.frame(matrix(nrow = num_SAEM, ncol = 3))
  for (i in 1:num_SAEM) {
    CI_SAEM[i, ] = (estimates_SAEM[i, ] - 1.96 * SE_SAEM[i, ] <= true_value & estimates_SAEM[i, ] + 1.96 * SE_SAEM[i, ] >= true_value)
  }
  # Find %rMSE, bias, and %bias
  # %rMSE = sqrt([\sum_{i=1}^num_SAEM (beta(i) - beta)^2]/num_SAEM)/beta
  # bias = [\sum_{i=1}^num_SAEM (beta(i) - beta)]/num_SAEM
  # %bias = bias/beta
  rMSE_SAEM <- bias_SAEM <- bias_perc_SAEM <- vector("double", length(all_parameters))
  for (i in seq_along(all_parameters)) {
    rMSE_SAEM[i] = sqrt(sum((estimates_SAEM_all[, i] - all_parameters[i]) ^ 2) / num_SAEM) / all_parameters[i]
    bias_SAEM[i] = sum(estimates_SAEM_all[, i] - all_parameters[i]) / num_SAEM
    bias_perc_SAEM[i] = bias_SAEM[i] / all_parameters[i]
  }
  Coverage_SAEM <- sapply(CI_SAEM, mean)
  result_SAEM <- cbind(true_value = all_parameters, 
                       estimates = sapply(estimates_SAEM_all, mean),
                       se_m = c(sapply(SE_SAEM, mean), rep(NA, 4)), # model SE
                       se_s = sapply(estimates_SAEM_all, sd), # simulation SE
                       rMSE_SAEM, # %rMSE
                       bias_SAEM, # bias
                       bias_perc_SAEM, # %bias
                       CP = c(Coverage_SAEM, rep(NA, 4)), # 95% CI coverage rate
                       NC = c(num_sim - num_SAEM, rep(NA, 6)), # number of non-convergence
                       time = c(sum(time_SAEM), rep(NA, 6))) # total computation time
  rownames(result_SAEM) = c("p1", "lambda1", "p2", "D11", "D22", "D33", "sigma")
  
  # Results table for NLME
  estimates_NLME_all <- cbind(estimates_NLME, D_NLME, sigma_NLME)
  # In case of any non-convergence, number of valid simulation iteration is
  num_NLME <- nrow(estimates_NLME)
  # CI_NLME is a data frame containing boolean, TRUE is the estimate is within 95% CI, FALSE otherwise.
  CI_NLME <- data.frame(matrix(nrow = num_NLME, ncol = 3))
  for (i in 1:num_NLME) {
    CI_NLME[i, ] = (estimates_NLME[i, ] - 1.96 * SE_NLME[i, ] <= true_value & estimates_NLME[i, ] + 1.96 * SE_NLME[i, ] >= true_value)
  }
  # Find %rMSE, bias, and %bias
  # %rMSE = sqrt([\sum_{i=1}^num_NLME (beta(i) - beta)^2]/num_NLME)/beta
  # bias = [\sum_{i=1}^num_NLME (beta(i) - beta)]/num_NLME
  # %bias = bias/beta
  rMSE_NLME <- bias_NLME <- bias_perc_NLME <- vector("double", length(all_parameters))
  for (i in seq_along(all_parameters)) {
    rMSE_NLME[i] = sqrt(sum((estimates_NLME_all[, i] - all_parameters[i]) ^ 2) / num_NLME) / all_parameters[i]
    bias_NLME[i] = sum(estimates_NLME_all[, i] - all_parameters[i]) / num_NLME
    bias_perc_NLME[i] = bias_NLME[i] / all_parameters[i]
  }
  Coverage_NLME <- sapply(CI_NLME, mean)
  result_NLME <- cbind(all_parameters,
                       estimates = sapply(estimates_NLME_all, mean),
                       se_m = c(sapply(SE_NLME, mean), rep(NA, 4)), # model SE
                       se_s = sapply(estimates_NLME_all, sd),  # simulation SE
                       rMSE_NLME, # %rMSE
                       bias_NLME, # bias
                       bias_perc_NLME, # %bias
                       CP = c(Coverage_NLME, rep(NA, 4)), # 95% CI coverage rate
                       NC = c(num_sim - num_NLME, rep(NA, 6)), # number of non-convergence
                       time = c(sum(time_NLME), rep(NA, 6))) # total computation time
  rownames(result_NLME) = c("p1", "lambda1", "p2", "D11", "D22", "D33", "sigma")
  
  # Results table for LME4
  estimates_LME4_all <- cbind(estimates_LME4, D_LME4, sigma_LME4)
  # In case of any non-convergence, number of valid simulation iteration is
  num_LME4 <- nrow(estimates_LME4)
  # CI_LME4 is a data frame containing boolean, TRUE is the estimate is within 95% CI, FALSE otherwise.
  CI_LME4 <- data.frame(matrix(nrow = num_LME4, ncol = 3))
  for (i in 1:num_LME4) {
    CI_LME4[i, ] = (estimates_LME4[i, ] - 1.96 * SE_LME4[i, ] <= true_value & estimates_LME4[i, ] + 1.96 * SE_LME4[i, ] >= true_value)
  }
  # Find %rMSE, bias, and %bias
  # %rMSE = sqrt([\sum_{i=1}^num_LME4 (beta(i) - beta)^2]/num_LME4)/beta
  # bias = [\sum_{i=1}^num_LME4 (beta(i) - beta)]/num_LME4
  # %bias = bias/beta
  rMSE_LME4 <- bias_LME4 <- bias_perc_LME4 <- vector("double", length(all_parameters))
  for (i in seq_along(all_parameters)) {
    rMSE_LME4[i] = sqrt(sum((estimates_LME4_all[, i] - all_parameters[i]) ^ 2) / num_LME4) / all_parameters[i]
    bias_LME4[i] = sum(estimates_LME4_all[, i] - all_parameters[i]) / num_LME4
    bias_perc_LME4[i] = bias_LME4[i] / all_parameters[i]
  }
  Coverage_LME4 <- sapply(CI_LME4, mean)
  result_LME4 <- cbind(all_parameters,
                       estimates = sapply(estimates_LME4_all, mean),
                       se_m = c(sapply(SE_LME4, mean), rep(NA, 4)), # model SE
                       se_s = sapply(estimates_LME4_all, sd), # simulation SE
                       rMSE_LME4, # %rMSE
                       bias_LME4, # bias
                       bias_perc_LME4, # %bias
                       CP = c(Coverage_LME4, rep(NA, 4)), # 95% CI coverage rate
                       NC = c(num_sim - num_LME4, rep(NA, 6)), # number of non-convergence
                       time = c(sum(time_LME4), rep(NA, 6)) # total computation time
  )
  rownames(result_LME4) = c("p1", "lambda1", "p2", "D11", "D22", "D33", "sigma")
  
  # Latex Table
  table_names <- c("Method", "Time (s)", "NC", "Parameter", "True Value", "Estimate", 
                   "SE_M", "SE_S", "Bias (%)", "rMSE (%)", "Coverage")
  table <- data.frame(matrix(nrow = 21, ncol = length(table_names)))
  colnames(table) <- table_names
  table$Parameter <- rep(c("$P_1$", "$\\lambda_1$", "$P_2$", "$D_{11}$", "$D_{22}$", "$D_{33}$", "$\\sigma_1$"), 3)
  table$`True Value` <- rep(all_parameters, 3)
  table$Method <- c("SAEM", rep(NA, 6), "NLME", rep(NA, 6), "LME4", rep(NA, 6))
  table$Estimate <- c(result_SAEM[, 2], result_NLME[, 2], result_LME4[, 2])
  table$SE_M <- c(result_SAEM[, 3], result_NLME[, 3], result_LME4[, 3])
  table$SE_S <- c(result_SAEM[, 4], result_NLME[, 4], result_LME4[, 4])
  table$`Bias (%)` <- c(result_SAEM[, 7] * 100, result_NLME[, 7] * 100, result_LME4[, 7] *
                          100)
  table$`rMSE (%)` <- c(result_SAEM[, 5] * 100, result_NLME[, 5] * 100, result_LME4[, 5] *
                          100)
  table$Coverage <-  c(result_SAEM[, 8] * 100, result_NLME[, 8] * 100, result_LME4[, 8] *
                         100)
  table$`Time (s)` <- c(result_SAEM[, 10], result_NLME[, 10], result_LME4[, 10])
  table$NC <- c(result_SAEM[, 9], result_NLME[, 9], result_LME4[, 9])
  write.csv(table, file = paste0("Simulation/", folder.name, "/table.csv"), row.names = FALSE, na = "")
  
  latex_table <- xtable(table, type = "latex", align = c("cccccccccccc"))
  digits(latex_table) <- c(0, 2, 1, 0, 1, 1, 2, 2, 2, 2, 2, 2)
  print(latex_table,
        file = paste0("Simulation/", folder.name, "/table.tex"),
        include.rownames = FALSE,
        sanitize.text.function = function(x) {x},
        hline.after = c(-1, 0, 21))
}