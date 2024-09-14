##=====================Decay Model=======================
# Decay model: y_ij = log10(exp(p_1i - lambda_1i * t_ij) + exp(p_2i)) + e_ij
# Random effects: p_1i = p1 + b_1i, p_2i = p1 + b_2i, lambda_1i = lambda1 + b_3i,
# where e_ij ~ N(0, sigma1^2) and b_i ~ N(0, D) with diagonal elements being D11, D22, and D33.

##=====================Setting======================
SimFunDecay <- function(p1_t = 17, lambda1_t = 4, p2_t = 3, D = diag(c(2, 0.1, 0.1)),
                        sigma = 0.1, time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 5, num_patient = 20) {
  # Find which parameter has a random effect
  D_diag <<- which(diag(D) != 0)
  D_name <- ""
  for (i in 1:length(D_diag)){
    D_name = paste(D_name, "_D", D_diag[i], D_diag[i], "_", D[D_diag[i], D_diag[i]], sep = "")
  }
  
  # Create a folder to store simulation results
  folder.name = paste("Decay_n", num_patient, "_ni", length(time), "_p1_", p1_t, "_lambda1_", lambda1_t,
                      "_p2_", p2_t, "_sigma", sigma, D_name, "_numsim", num_sim, sep="")
  dir.create(file.path("Simulation/", folder.name), showWarnings = FALSE)
  
  # Parameter estimates
  estimates_SAEM <- estimates_NLME <- estimates_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(estimates_SAEM) <- colnames(estimates_NLME) <- colnames(estimates_LME4) <- c("p1", "lambda1", "p2")
  # SE of the parameter estimates
  SE_SAEM <- SE_NLME <- SE_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(SE_SAEM) <- colnames(SE_NLME) <- colnames(SE_LME4) <- c("p1", "lambda1", "p2")
  # Variance of random effects
  D_SAEM <- D_NLME <- D_LME4 <- data.frame(matrix(nrow = num_sim, ncol = length(D_diag)))
  # Residual SD (sigma_1)
  sigma_SAEM <- sigma_NLME <- sigma_LME4 <- c(rep(NA, num_sim))
  # Computation time
  time_SAEM <- time_NLME <- time_LME4 <- c(rep(NA, num_sim))

  D_names <- c()
  for (i in 1:length(D_diag)){
    D_names[i] = paste("D", D_diag[i], D_diag[i], sep = "")
  }
  colnames(D_SAEM) <- colnames(D_NLME) <- colnames(D_LME4) <- D_names
  
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
      mutate(viral = log10(exp(p1_t + b1_sim - (lambda1_t + b3_sim) * day) + exp(p2_t + b2_sim)) + e) %>%
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
    # Find which parameter has a random effect
    D_rand <- diag(D) != 0
    setIndividualParameterVariability(p1 = D_rand[1], p2 = D_rand[2], lambda1 = D_rand[3])
    
    # Set the population parameter information
    setPopulationParameterInformation(p1_pop = list(initialValue = p1_t),
                                      p2_pop = list(initialValue = p2_t),
                                      lambda1_pop = list(initialValue = lambda1_t))
    setIndividualParameterDistribution(p1 = "normal", lambda1 = "normal", p2 = "normal")
    setScenario(scenario)
    # Run the estimation
    runScenario()
    
    # Store the estimates in table
    estimates_SAEM[i, ] <- getEstimatedPopulationParameters()[c(1, 3, 2)]
    SE_SAEM[i, ] <- getEstimatedStandardErrors()$stochasticApproximation[c(1, 3, 2), 2]
    D_SAEM[i, ] <- getEstimatedPopulationParameters()[4: (length(D_diag)+3)] 
    sigma_SAEM[i] <- tail(getEstimatedPopulationParameters(), 1)
    end.time <- Sys.time()
    time_SAEM[i] <- end.time - start.time
    
    # ================ Fit viral decay model using NLME =================
    data_nlme <- groupedData(viral ~ day | PATIENT, data = data_sim)
  
    start.time <- Sys.time()
    randef <<- c("p1", "p2", "lambda1")[D_rand]
    randef_name <- paste(paste0(randef), collapse = " + ")
    
    tryCatch({
      model_nlme <- nlme(
        viral ~ log10(exp(p1 - lambda1 * day) + exp(p2)),
        fixed = p1 + lambda1 + p2 ~ 1,
        random = as.formula(paste0(randef_name, "~", "1")),
        data = data_nlme,
        start = c(p1 = p1_t, b1 = lambda1_t, p2 = p2_t)
      )
      
      estimates_NLME[i, ] <- summary(model_nlme)$tTable[, 1]
      SE_NLME[i, ] <- summary(model_nlme)$tTable[, 2]
      D_NLME[i, ] <- as.numeric(VarCorr(model_nlme)[1:length(D_diag), 1])
      sigma_NLME[i] <- as.numeric(VarCorr(model_nlme)[length(D_diag)+1, 2])
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
    # Dynamically construct the random effects part
    model_lme4 <- nlmer(
      formula = as.formula(paste("viral ~ nfun(day, p1, lambda1, p2) ~", paste(paste0("(", randef, " | PATIENT)"), collapse = " + "))),
      data_sim,
      start = c(p1 = p1_t, lambda1 = lambda1_t, p2 = p2_t)
    )
    
    estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
    SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
    D_LME4[i, ] <- as.numeric(as.data.frame(VarCorr(model_lme4))[1:length(D_diag), 4])
    sigma_LME4[i] <- as.numeric(as.data.frame(VarCorr(model_lme4))[length(D_diag)+1, 5])
    }, error = function(e) {
      estimates_LME4[i, ] <<- SE_LME4[i, ] <<- D_LME4[i, ] <<- sigma_LME4[i] <<- NA
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
                D_SAEM = D_SAEM, D_NLME = D_NLME, D_LME4 = D_LME4)
  files_sigma <- list(sigma_SAEM = sigma_SAEM, sigma_NLME = sigma_NLME, sigma_LME4 = sigma_LME4)
  files_time <- list(time_SAEM = time_SAEM, time_NLME = time_NLME, time_LME4 = time_LME4)
  
  # Remove NA's
  files_est <- lapply(files_est, function(df) df[!apply(is.na(df), 1, all), ])
  files_sigma <- lapply(files_sigma, na.omit)
  files <- c(files_est, files_sigma, files_time)
  # Save simulation results
  for (i in seq_along(files)) {
    write.table(files[[i]], file = paste0("Simulation/", folder.name, "/Results/", names(files)[i], ".csv"), row.names = FALSE)
  }
  
  #======================Calculate MSE and bias===================
  true_value <- c(p1_t, lambda1_t, p2_t)
  all_parameters <- c(true_value, diag(D)[D_diag], sigma)

  # Results table for SAEM
  result_SAEM <- CreateResultTable(num_sim, TrueValue = true_value, TrueCov = D, 
                                   TrueSigma = sigma, EstTable = files$estimates_SAEM, SETable = files$SE_SAEM, 
                                   CovTable = files$D_SAEM, SigmaTable = files$sigma_SAEM, TimeTable = files$time_SAEM)
  # Results table for NLME
  result_NLME <- CreateResultTable(num_sim, true_value, D, sigma, files$estimates_NLME, files$SE_NLME, files$D_NLME, files$sigma_NLME, files$time_NLME)
  
  # Results table for LME4
  result_LME4 <- CreateResultTable(num_sim, true_value, D, sigma, files$estimates_LME4, files$SE_LME4, files$D_LME4, files$sigma_LME4, files$time_LME4)

  resultTable <- rbind(result_SAEM, result_NLME, result_LME4)

  # Latex Table
  D_latex <- c()
  for(i in 1:length(D_diag)){
    D_latex <- c(D_latex, paste("$", "D_{", D_diag[i], D_diag[i], "}$", sep = ""))
  }
  table <- CreateLatexTable(resultTable, 
                            ParameterName = c("$P_1$", "$\\lambda_1$", "$P_2$"), 
                            CovName = D_latex, 
                            SigmaName = "$\\sigma_1$")
  
  write.csv(table, file = paste0("Simulation/", folder.name, "/table.csv"), row.names = FALSE, na = "")
  
  latex_table <- xtable(table, type = "latex", align = c("cccccccccccc"))
  digits(latex_table) <- c(0, 2, 1, 0, 1, 1, 2, 2, 2, 2, 2, 2)
  print(latex_table,
        file = paste0("Simulation/", folder.name, "/table.tex"),
        include.rownames = FALSE,
        sanitize.text.function = function(x) {x},
        hline.after = c(-1, 0, length(all_parameters) * 3))
}