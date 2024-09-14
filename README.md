# Simulation Studies for NLME and GLMM Models

This repository contains simulation studies and model fitting scripts for **Nonlinear Mixed Effects Models (NLME)** and **Generalized Linear Mixed Models (GLMM)**. The simulations are primarily used for assessing model performance in various contexts, including population models, viral decay models, and handling of random effects.

## Overview

The repository includes:

- **Viral Load Decay Simulation**: Simulates viral load data based on a decay model before ART interruption, including random effects for patient-specific parameters.
- **Model Fitting**: Implements model fitting using different methods like SAEM (Stochastic Approximation Expectation-Maximization), NLME (Nonlinear Mixed Effects), and LME4 (Linear and Nonlinear Mixed Models).
- **Performance Metrics**: Measures and compares estimation performance (bias, RMSE, SE) across methods and evaluates computation time, convergence, and confidence interval coverage.

### Key Functions

1. **SimFunDecay()**: 
   - Simulates viral load decay data using patient-specific parameters.
   - Fits NLME models using various techniques (SAEM, NLME, LME4).
   - Stores and compares estimates, standard errors, and variance components.

2. **CreateResultTable()**: 
   - Summarizes model performance for parameters, random effects, and residual variance.
   - Computes bias, RMSE, and confidence interval coverage for parameter estimates.

3. **CreateLatexTable()**:
   - Generates a LaTeX-ready table summarizing the results of the simulation study.

## Simulation Settings

- **Model**: Viral load decay model: \( y_{ij} = \log_{10} (\exp(p_{1i} - \lambda_{1i} t_{ij}) + \exp(p_{2i})) + e_{ij} \)
- **Parameters**:
  - \( p_{1i}, \lambda_{1i}, p_{2i} \): Patient-specific parameters with random effects
  - \( D \): Covariance matrix for random effects
  - \( \sigma \): Residual standard deviation
- **Methods**: 
  - SAEM for maximum likelihood estimation
  - NLME using grouped data and random effects
  - LME4 for nonlinear mixed effects modeling

## Prerequisites

- R packages: `nlme`, `lme4`, `MASS`, `ggplot2`, `dplyr`
- Monolix Suite for SAEM estimation (via `lixoftConnector`)

## Usage

Clone the repository:

```bash
git clone https://github.com/UBCSG/Simulation_NLME_GLMM.git
```

Run simulations with:

```r
SimFunDecay(num_sim = 5, num_patient = 20)
```

The results will be saved in the `Simulation/` directory, with individual results for each simulation method and a final LaTeX table summarizing the results.

## Contributing

Feel free to submit issues or pull requests for improvements, bug fixes, or additional model types.

## License

This project is licensed under the MIT License.

---

You can customize this further based on your repository's content and focus.
