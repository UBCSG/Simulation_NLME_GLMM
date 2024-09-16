##=====================Rebound Model=======================
# Rebound model: y_ij = beta_1i * t_ij / (t_ij  + exp( beta_2i - beta_3i * t_ij)) + beta_4i / (1 + exp( beta_5i * t)) + e_ij
# Random effects: beta_ki = betak + tau_ki, k = 1, 2, 3, 4, 5
# where e_ij ~ N(0, sigma2^2) and tau_i ~ N(0, B) with diagonal elements being B11, B22, B33, and B44.

##=====================Setting======================
SimFunRebound <- function(beta1_t = 4, beta2_t = 5.7, beta3_t = 2.1, beta4_t = 1.9, beta5_t = 0.4,
                          B = diag(c(0.2, 0, 0.1, 0, 0)), sigma = 0.1, 
                          time = c(0.5, 2.9, 4.8, 7, 10.2), num_sim = 500, num_patient = 50) {
  
  # Find which parameter has a random effect
  B_diag <<- which(diag(B) != 0)
  B_name <- ""
  for (i in 1:length(B_diag)){
    B_name = paste(B_name, "_B", B_diag[i], B_diag[i], "_", B[B_diag[i], B_diag[i]], sep = "")
  }

  # Create a folder to store simulation results
  folder.name = paste("Rebound_n", num_patient, "_ni", length(time), "_beta1_", beta1_t, 
                                    "_beta2_", beta2_t, "_beta3_", beta3_t, "_beta4_", beta4_t, "_beta5_", beta5_t, 
                                    "_sigma", sigma,  B_name,"_numsim", num_sim, sep = "")
  dir.create(file.path("Simulation/", folder.name), showWarnings = FALSE)
  
  # Parameter estimates
  estimates_SAEM <- estimates_NLME <- estimates_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 5))
  colnames(estimates_SAEM) <- colnames(estimates_NLME) <- colnames(estimates_LME4) <- c("beta1", "beta2", "beta3", "beta4", "beta5")
  # SE of the parameter estimates
  SE_SAEM <- SE_NLME <- SE_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 5))
  colnames(SE_SAEM) <- colnames(SE_NLME) <- colnames(SE_LME4) <- c("beta1", "beta2", "beta3", "beta4", "beta5")
  # Variance of random effects
  B_SAEM <- B_NLME <- B_LME4 <- data.frame(matrix(nrow = num_sim, ncol = length(B_diag)))
  # Residual SD (sigma2)
  sigma_SAEM <- sigma_NLME <- sigma_LME4 <- c(rep(NA, num_sim))
  # Computation time
  time_SAEM <- time_NLME <- time_LME4 <- c(rep(NA, num_sim))
  
  B_names <- c()
  for (i in 1:length(B_diag)){
    B_names[i] = paste("B", B_diag[i], B_diag[i], sep = "")
  }
  colnames(B_SAEM) <- colnames(B_NLME) <- colnames(B_LME4) <- B_names
  
  for (i in 1:num_sim) {
    # Print the iteration number
    print(glue::glue("i=", i))
    
    # Simulate dataset for viral load after ART interruption
    # Create a vector 'PATIENT' where each patient ID is repeated for each time point
    PATIENT <- rep(1:num_patient, each = length(time))
    # Create a vector 'day' that repeats the time points for each patient
    day <- rep(time, num_patient)
    # Combine the 'PATIENT' and 'day' vectors into a matrix 'data'
    data <- cbind(PATIENT, day)
    
    # Simulate b_i and e_ij
    bi_sim <- mvrnorm(num_patient, c(0, 0, 0, 0, 0), B)
    colnames(bi_sim) = c("b1_sim", "b2_sim", "b3_sim", "b4_sim", "b5_sim")
    bi_sim <- cbind(PATIENT = c(1:num_patient), bi_sim)
    data <- merge(data, bi_sim, by = "PATIENT")
    data <- data %>%
      mutate(e = rnorm(nrow(data), 0, sigma))
    
    # Simulate viral load based on y_ij = beta_1i * t_ij / (t_ij  + exp( beta_2i - beta_3i * t_ij)) + beta_4i / (1 + exp( beta_5i * t)) + e_ij
    data_sim <- data %>%
      mutate(viral = (beta1_t + b1_sim) * day / (day  + exp(beta2_t + b2_sim - (beta3_t + b3_sim) * day)) + (beta4_t + b4_sim) / (1 + exp((beta5_t + b5_sim) * day)) + e) %>%
      dplyr::select(-b1_sim, -b2_sim, -b3_sim, -b4_sim, -b5_sim, -e)
   
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
    
    # ================ Fit viral rebound model using SAEM =================
    start.time <- Sys.time()
    data = list(
      dataFile = paste0("Simulation/", folder.name, "/data.txt"),
      headerTypes = c("id", "time", "observation"),
      observationTypes = list(viral = "continuous")
    )
    modelFile = paste0('Model/model_rebound.txt')
    newProject(data = data, modelFile = modelFile)
    # Set the error model type to be used for the observation model
    setErrorModel(viral = "constant") # obs = pred + a*err, error ~ N(0, 1)
    # Set observation model distribution
    setObservationDistribution(viral = "normal") 
    # Set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation = T,
                       logLikelihoodEstimation = T)
    scenario$linearization = FALSE
    # Find which parameter has a random effect
    B_rand <- diag(B) != 0
    setIndividualParameterVariability(beta1 = B_rand[1], beta2 = B_rand[2], 
                                      beta3 = B_rand[3], beta4 = B_rand[4], beta5 = B_rand[5])

    # Set the population parameter information
    setPopulationParameterInformation(beta1_pop = list(initialValue = beta1_t),
                                      beta2_pop = list(initialValue = beta2_t),
                                      beta3_pop = list(initialValue = beta3_t),
                                      beta4_pop = list(initialValue = beta4_t),
                                      beta5_pop = list(initialValue = beta5_t))
    setIndividualParameterDistribution(beta1 = "normal", beta2 = "normal", 
                                       beta3 = "normal", beta4 = "normal", beta5 = "normal")
    setScenario(scenario)
    # Run the estimation
    runScenario()
    
    # Store the estimates in table
    estimates_SAEM[i, 1:5] <- getEstimatedPopulationParameters()[1:5]
    SE_SAEM[i, ] <- getEstimatedStandardErrors()$stochasticApproximation[1:5, 2]
    sigma_SAEM[i] <- tail(getEstimatedPopulationParameters(), 1)
    B_SAEM[i, ] <- (getEstimatedPopulationParameters()[6:(5+length(B_diag))])^2
    end.time <- Sys.time()
    time_SAEM[i] <- end.time - start.time
    
    # ================ Fit viral rebound model using NLME =================
    data_nlme <- groupedData(viral ~ day | PATIENT, data = data_sim)
    
    start.time <- Sys.time()
    beta_name <- paste(paste0("beta", B_diag, sep = ""), collapse = " + ")
    
    tryCatch({
      model_nlme <- nlme(
        viral ~ beta1 * day / (day  + exp(beta2 - beta3 * day)) + beta4 / (1 + exp(beta5 * day)),
        fixed = beta1 + beta2 + beta3 + beta4 + beta5 ~ 1,
        random = as.formula(paste0(beta_name, "~", "1")) ,
        data = data_nlme,
        start = c(beta1 = beta1_t, beta2 = beta2_t, beta3 = beta3_t, beta4 = beta4_t, beta5 = beta5_t)
      )
      
      estimates_NLME[i, ] <- summary(model_nlme)$tTable[, 1]
      SE_NLME[i, ] <- summary(model_nlme)$tTable[, 2]
      B_NLME[i, ] <- as.numeric(VarCorr(model_nlme)[1:length(B_diag), 1])
      sigma_NLME[i] <- as.numeric(VarCorr(model_nlme)[length(B_diag) + 1, 2])
    }, error = function(e) {
      estimates_NLME[i, ] <<- SE_NLME[i, ] <<- B_NLME[i, ] <<- sigma_NLME[i]  <<- NA
    })
    end.time <- Sys.time()
    time_NLME[i] <- end.time - start.time
    
    
    # ================ Fit viral rebound model using LME4 =================
    start.time <- Sys.time()
    
    nform <- ~ (beta1 * t / (t  + exp(beta2 - beta3 * t)) + beta4 / (1 + exp(beta5 * t)))
    nfun <- deriv(
      nform,
      namevec = c("beta1", "beta2", "beta3", "beta4", "beta5"),
      function.arg = c("t", "beta1", "beta2", "beta3", "beta4", "beta5")
    )
    
    tryCatch({
      # Dynamically construct the random effects part
      model_lme4 <- nlmer(
        formula = as.formula(paste("viral ~ nfun(day, beta1, beta2, beta3, beta4, beta5) ~",  paste(paste0("(beta", B_diag, " | PATIENT)"), collapse = " + "))),
        data = data_nlme,
        start = c(beta1 = beta1_t, beta2 = beta2_t, beta3 = beta3_t, beta4 = beta4_t, beta5 = beta5_t)
      )
      
      estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
      SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
      B_LME4[i, ] <- as.numeric(as.data.frame(VarCorr(model_lme4))[1:length(B_diag), 4])
      sigma_LME4[i] <- as.numeric(as.data.frame(VarCorr(model_lme4))[length(B_diag)+1, 5])
    }, error = function(e) {
      estimates_LME4[i, ] <<- SE_LME4[i, ] <<- B_LME4[i, ] <<- sigma_LME4[i] <<- NA 
    })
    
    end.time <- Sys.time()
    time_LME4[i] <- end.time - start.time
  }
  
  # Create a folder to store the simulation results
  dir.create(file.path("Simulation/", folder.name, "Results"), showWarnings = FALSE)
  SE_SAEM[SE_SAEM == "NaN"] <- NA
  SE_SAEM <- mutate_all(SE_SAEM, function(x) as.numeric(as.character(x)))
  files_est <- list(estimates_SAEM = estimates_SAEM, estimates_NLME = estimates_NLME, estimates_LME4 = estimates_LME4,
                SE_SAEM = SE_SAEM, SE_NLME = SE_NLME, SE_LME4 = SE_LME4,
                B_SAEM = B_SAEM, B_NLME = B_NLME, B_LME4 = B_LME4)
  files_sigma <- list(sigma_SAEM = sigma_SAEM, sigma_NLME = sigma_NLME, sigma_LME4 = sigma_LME4)
  files_time <- list(time_SAEM = time_SAEM, time_NLME = time_NLME, time_LME4 = time_LME4)
              
  # Remove NA's
  files_est <- lapply(files_est, function(df) df[!apply(is.na(df), 1, all), ])
  files_sigma <- lapply(files_sigma, na.omit)
  files <- c(files_est, files_sigma, files_time)
  # Save simulation results
  for (i in seq_along(files)) {
    write.csv(files[[i]], file = paste0("Simulation/", folder.name, "/Results/", names(files)[i], ".csv"), row.names = FALSE)
  }
  
  #======================Calculate MSE and bias===================
  true_value <- c(beta1_t, beta2_t, beta3_t, beta4_t, beta5_t)
  all_parameters <- c(true_value, diag(B)[B_diag], sigma)
  
  # Results table for SAEM
  result_SAEM <- CreateResultTable(num_sim, true_value, B, sigma, files$estimates_SAEM, files$SE_SAEM, files$B_SAEM, files$sigma_SAEM, files$time_SAEM)

  # Results table for NLME
  result_NLME <- CreateResultTable(num_sim, true_value, B, sigma, files$estimates_NLME, files$SE_NLME, files$B_NLME, files$sigma_NLME, files$time_NLME)

  # Results table for LME4
  result_LME4 <- CreateResultTable(num_sim, true_value, B, sigma, files$estimates_LME4, files$SE_LME4, files$B_LME4, files$sigma_LME4, files$time_LME4)

  resultTable <- rbind(result_SAEM, result_NLME, result_LME4)
  
  # Latex Table
  B_latex <- c()
  for(i in 1:length(B_diag)){
    B_latex <- c(B_latex, paste("$", "B_{", B_diag[i], B_diag[i], "}$", sep = ""))
  }
  table <- CreateLatexTable(resultTable, 
                            ParameterName = c("$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$", "$\\beta_5$"), 
                            CovName = B_latex, 
                            SigmaName = "$\\sigma_2$")

  write.csv(table, file = paste0("Simulation/", folder.name, "/table.csv"), row.names = FALSE, na = "")
  
  latex_table <- xtable(table, type = "latex", align = c("cccccccccccc"))
  digits(latex_table) <- c(0, 2, 1, 0, 1, 1, 2, 2, 2, 2, 2, 2)

  print(latex_table,
        file = paste0("Simulation/", folder.name, "/table.tex"),
        include.rownames = FALSE,
        sanitize.text.function = function(x) {x},
        hline.after = c(-1, 0, 3 * (length(all_parameters))))
}
