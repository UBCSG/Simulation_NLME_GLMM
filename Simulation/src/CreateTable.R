CreateResultTable <- function(num_sim, TrueValue, TrueCov, Cov_names, TrueSigma, EstTable, SETable, CovTable, SigmaTable, TimeTable) {
  randef <- which(diag(TrueCov) != 0)
  
  all_parameters <- c(TrueValue, diag(TrueCov)[randef], TrueSigma)
  estimates_all <- cbind(EstTable, CovTable, SigmaTable)
  # In case of any non-convergence, number of valid simulation iteration is
  num_iter <- nrow(EstTable)
  # CI_SAEM is a data frame containing boolean, TRUE is the estimate is within 95% CI, FALSE otherwise.
  CI_table <- data.frame(matrix(nrow = num_iter, ncol = ncol(EstTable)))
  for (i in 1:num_iter){
    CI_table[i, ] = (EstTable[i, ] - 1.96 * SETable[i, ] <= TrueValue & EstTable[i, ] + 1.96 * SETable[i, ] >= TrueValue)
  }
  # Find %rMSE, bias, and %bias
  # %rMSE = sqrt([\sum_{i=1}^num_iter (beta(i) - beta)^2]/num_iter)/beta
  # bias = [\sum_{i=1}^num_iter (beta(i) - beta)]/num_iter
  # %bias = bias/beta
  rMSE <- bias <- bias_perc <- vector("double", length(all_parameters))
  for (i in seq_along(all_parameters)){
    rMSE[i] = sqrt(sum((estimates_all[, i] - all_parameters[i]) ^ 2) / num_iter) / all_parameters[i]
    bias[i] = sum(estimates_all[, i] - all_parameters[i]) / num_iter
    bias_perc[i] = bias[i] / all_parameters[i]
  }
  Coverage <- sapply(CI_table, function(x) mean(x, na.rm = TRUE))
  result <- cbind(true_value = all_parameters, 
                  estimates = sapply(estimates_all, mean), 
                  se_m = c(sapply(SETable, function(x) mean(x, na.rm = TRUE)), rep(NA, length(randef) + 1)), # model SE
                  se_s = sapply(estimates_all, sd), # simulation SE
                  rMSE, # %rMSE
                  bias, # bias
                  bias_perc, # %bias
                  CP = c(Coverage, rep(NA, length(randef) + 1)), # 95% CI coverage rate
                  NC = c(num_sim - num_iter, rep(NA, length(randef) + 5)), # number of non-convergence
                  time = c(sum(TimeTable), rep(NA, length(randef) + 5))) # total computation time
  rownames(result) = c("beta1", "beta2", "beta3", "beta4", "beta5", Cov_names, "sigma")
  as.data.frame(result)
}

CreateLatexTable <- function(resultTable, ParameterName, CovName, SigmaName) {
  table_names <- c("Method", "Time (s)", "NC", "Parameter", "True Value", "Estimate", 
                   "SE_M", "SE_S", "Bias (%)", "rMSE (%)", "Coverage")
  
  table <- data.frame(matrix(nrow = 3*(length(CovName)+6), ncol = length(table_names)))
  colnames(table) <- table_names

  table$Parameter <- rep(c(ParameterName, CovName, SigmaName), 3)
  table$`True Value` <- resultTable$true_value
  table$Method <- c("SAEM", rep(NA, length(CovName) + 5), "NLME", rep(NA, length(CovName) + 5), "LME4", rep(NA, length(CovName) + 5))
  table$Estimate <- resultTable$estimates
  table$SE_M <- resultTable$se_m
  table$SE_S <- resultTable$se_s
  table$`Bias (%)` <- resultTable$bias_perc * 100
  table$`rMSE (%)` <- resultTable$rMSE * 100
  table$Coverage <-  resultTable$CP
  table$`Time (s)` <- resultTable$time
  table$NC <- resultTable$NC
  table
}
