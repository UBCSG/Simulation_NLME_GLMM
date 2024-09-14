##=====================GLMM=======================
# GLMM: logit(y_ij) = (alpha_0 + a_{0i}) + (alpha_1 + a_{1i}) * t_{ij} + (alpha_2 + a_{2i}) * x_{ij}, 
# where x_{ij} ~ N(0, 1) is a covariate, and a_i ~ N(0, A) with diagonal elements being A11, A22, and A33.

##=====================Setting======================
SimFunGLMM <- function(alpha0_t = 8, alpha1_t = -2, alpha2_t = 2, A = diag(c(2, 0.5, 0)),
                          time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 5, num_patient = 50) {
  # Find which parameter has a random effect
  A_diag <<- which(diag(A) != 0)
  A_name <- ""
  for (i in 1:length(A_diag)){
    A_name = paste(A_name, "_A", A_diag[i], A_diag[i], "_", A[A_diag[i], A_diag[i]], sep = "")
  }
  
  # Create a folder to store simulation results
  folder.name = paste("GLMM_n", num_patient, "_ni", length(time), "_alpha0_", alpha0_t, 
                      "_alpha1_", alpha1_t, "_alpha2_", alpha2_t, A_name, "_numsim", num_sim, sep = "")
  dir.create(file.path("Simulation/", folder.name), showWarnings = FALSE)
  
  # Parameter estimates
  estimates_SAEM <- estimates_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(estimates_SAEM) <- colnames(estimates_LME4) <- c("alpha0", "alpha1", "alpha2")
  # SE of the parameter estimates
  SE_SAEM <- SE_LME4 <- data.frame(matrix(nrow = num_sim, ncol = 3))
  colnames(SE_SAEM) <- colnames(SE_LME4) <- c("alpha0", "alpha1", "alpha2")
  # Variance of random effects
  A_SAEM <- A_LME4 <- data.frame(matrix(nrow = num_sim, ncol = length(A_diag)))
  # Computation time
  time_SAEM <- time_LME4 <- c(rep(NA, num_sim))
  
  A_names <- c()
  for (i in 1:length(A_diag)){
    A_names[i] = paste("A", A_diag[i], A_diag[i], sep = "")
  }
  colnames(A_SAEM) <- colnames(A_LME4) <- A_names
  
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
    
    # Simulate a_i
    ai_sim <- mvrnorm(num_patient, c(0, 0, 0), A)
    colnames(ai_sim) = c("a0_sim", "a1_sim", "a2_sim")
    ai_sim <- cbind(PATIENT = c(1:num_patient), ai_sim)
    data <- merge(data, ai_sim, by = "PATIENT")
    
    # Simulate response based on logit(y_ij) = (alpha_0 + a_{0i}) + (alpha_1 + a_{1i}) * t_{ij} + (alpha_2 + a_{2i}) * x_{ij}
    data_sim <- data %>%
      mutate(x = rnorm(nrow(data), 0, 1)) %>%
      mutate(xb = (alpha0_t + a0_sim) + (alpha1_t + a1_sim) * day + (alpha2_t + a2_sim) * x) %>%
      mutate(p = 1 / (1 + exp(-xb))) %>%
      mutate(y = rbinom(n = nrow(data), size = 1, p = p)) %>%
      dplyr::select(-a0_sim, -a1_sim, -a2_sim, -xb, -p)
    
    # Save the simulated data locally
    write.table(data_sim, paste0("Simulation/", folder.name, "/data.txt"), sep = "," , quote = FALSE, row.names = FALSE)
    
    # ================ Fit GLMM using SAEM =================
    start.time <- Sys.time()
    data = list(
      dataFile = paste0("Simulation/", folder.name, "/data.txt"),
      headerTypes = c("id", "time", "regressor", "observation"),
      observationTypes = list(y = "categorical")
    )
    modelFile = paste0('Model/model_GLMM.txt')
    newProject(data = data, modelFile = modelFile)
    
    # Set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation = T,
                       logLikelihoodEstimation = T)
    scenario$linearization = FALSE
    # Find which parameter has a random effect
    A_rand <- diag(A) != 0
    setIndividualParameterVariability(alpha0 = A_rand[1], alpha1 = A_rand[2], alpha2 = A_rand[3])
    setIndividualParameterDistribution(alpha0 = "normal", alpha1 = "normal", alpha2 = "normal")
    setPopulationParameterInformation(alpha0_pop = list(initialValue = alpha0_t),
                                      alpha1_pop = list(initialValue = alpha1_t),
                                      alpha2_pop = list(initialValue = alpha2_t))
    
    setScenario(scenario)
    # Run the estimation
    runScenario()
    
    # Store the estimates in table
    estimates_SAEM[i, ] <- getEstimatedPopulationParameters()[1:3]
    SE_SAEM[i, ] <- getEstimatedStandardErrors()$stochasticApproximation[1:3, 2]
    A_SAEM[i, ] <- getEstimatedPopulationParameters()[4:(3+length(A_diag))]
    end.time <- Sys.time()
    time_SAEM[i] <- end.time - start.time
    
    # ================ Fit GLMM using LME4 =================
    data_lme4 <- groupedData(y ~ day | PATIENT, data = data_sim)
    randef <<- c("1", "day", "x")[A_rand]
    ranef_names <<- paste(randef, collapse = "+")
    
    start.time <- Sys.time()
    tryCatch({
      model_lme4 <- glmer(
        formula = as.formula(paste("y ~ day + x + ",  paste("(", ranef_names, "|PATIENT)"))),
        data = data_lme4,
        family = binomial,
        start = list(fixef = c(alpha0 = alpha0_t, alpha1 = alpha1_t, alpha2 = alpha2_t))
      )
      conv <- summary(model_lme4)$optinfo$conv$lme4$code
      if (is.integer(conv)){
        estimates_LME4[i, ] <- SE_LME4[i, ] <- NA
      } else{
        estimates_LME4[i, ] <- summary(model_lme4)$coefficients[, 1]
        SE_LME4[i, ] <- summary(model_lme4)$coefficients[, 2]
        A_LME4[i, ] <- as.numeric(as.data.frame(VarCorr(model_lme4))[1:length(A_diag), 4])
      }
    }, error = function(e) {
      estimates_LME4[i, ] <<- SE_LME4[i, ] <<- A_LME4 <<- NA
    })
    end.time <- Sys.time()
    time_LME4[i] <- end.time - start.time
  }  
  
  # Create a folder to store the simulation results
  dir.create(file.path("Simulation/", folder.name, "Results"), showWarnings = FALSE)
  SE_SAEM[SE_SAEM == "NaN"] <- NA
  SE_SAEM <- mutate_all(SE_SAEM, function(x) as.numeric(as.character(x)))
  files_est <- list(estimates_SAEM = estimates_SAEM, estimates_LME4 = estimates_LME4,
                    SE_SAEM = SE_SAEM, SE_LME4 = SE_LME4,
                    A_SAEM = A_SAEM, A_LME4 = A_LME4)
  files_time <- list(time_SAEM = time_SAEM, time_LME4 = time_LME4)
  
  # Remove NA's
  files_est <- lapply(files_est, function(df) df[!apply(is.na(df), 1, all), ])
  files <- c(files_est, files_time)
  # Save simulation results
  for (i in seq_along(files)) {
    write.csv(files[[i]], file = paste0("Simulation/", folder.name, "/Results/", names(files)[i], ".csv"), row.names = FALSE)
  }
  
  #======================Calculate MSE and bias===================
  true_value <- c(alpha0_t, alpha1_t, alpha2_t)
  all_parameters <- c(true_value, diag(A)[A_diag])
  
  # Results table for SAEM
  result_SAEM <- CreateResultTable(num_sim, true_value, A, NULL, files$estimates_SAEM, files$SE_SAEM, files$A_SAEM, NULL, files$time_SAEM)

  # Results table for LME4
  result_LME4 <- CreateResultTable(num_sim, true_value, A, NULL, files$estimates_LME4, files$SE_LME4, files$A_LME4, NULL, files$time_LME4)
  
  resultTable <- rbind(result_SAEM, result_LME4)
  
  # Latex Table
  A_latex <- c()
  for(i in 1:length(A_diag)){
    A_latex <- c(A_latex, paste("$", "A_{", A_diag[i], A_diag[i], "}$", sep = ""))
  }
  table <- CreateLatexTable(resultTable, 
                            ParameterName = c("$\\alpha_0$", "$\\alpha_1$", "$\\alpha_2$"), 
                            CovName = A_latex, 
                            SigmaName = NULL)
  
  write.csv(table, file = paste0("Simulation/", folder.name, "/table.csv"), row.names = FALSE, na = "")
  
  latex_table <- xtable(table, type = "latex", align = c("cccccccccccc"))
  digits(latex_table) <- c(0, 2, 1, 0, 1, 1, 2, 2, 2, 2, 2, 2)
  
  print(latex_table,
        file = paste0("Simulation/", folder.name, "/table.tex"),
        include.rownames = FALSE,
        sanitize.text.function = function(x) {x},
        hline.after = c(-1, 0, 2 * (length(all_parameters))))
}


