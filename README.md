# Comparison of Computationally Efficient Approximate Methods for Nonlinear and Generalized Linear Mixed Effects Models

Generalized linear mixed models (GLMMs) and nonlinear mixed effects (NLME) models are popular in the analysis of longitudinal or clustered data. Statistical inference is typically based on likelihood methods. When the number of random effects in the models is large, the observed-data likelihood function involves high-dimensional and intractable integration, as these models are nonlinear in the (unobserved) random effects. "Exact" methods, such as Monte Carlo EM (MCEM) algorithms and numerical integration methods, can be computationally very intensive and may offer convergence issues. Computationally more efficient approximate methods, such as the stochastic approximation EM (SAEM) algorithm, linearization methods, or Laplace approximation methods, are therefore commonly used in practice. This repository contains the simulation code and functions used for evaluating the performance of these three methods for NLME and GLMM models. The SAEM, linearization, and Laplace approximation methods were implemented using the R packages `lixoftConnectors`, `nlme`, and `lme4`, respectively.

## Models Overview

- ** NLME model 1: a viral decay model during an antiretroviral treatment (Wu and Ding, 1999)**

\[
y_{ij} = \log_{10} \left( e^{P_{1i} - \lambda_{1i} t_{ij}} + e^{P_{2i}} \right) + e_{ij}, 
\]
where:
\[
P_{1i} = P_1 + b_{1i}, \quad P_{2i} = P_2 + b_{2i}, \quad \lambda_{1i} = \lambda_1 + b_{3i}, 
\]
\[
e_{ij} \sim N(0, \sigma_1^2), \quad \mathbf{b}_i \sim N(0, D).
\]

- **NLME model 2: a viral rebound model after treatment interruption (Wang et al., 2020)**

\[
y_{ij} = \dfrac{\beta_{1i} t_{ij}}{t  + \exp{(\beta_{2i} - \beta_{3i} t_{ij})}} + \dfrac{\beta_{4i}}{1 + \exp( \beta_{5i} t_{ij})} + e_{ij}, 
\]
where:
\[
\beta_{1i} = \beta_1 + \tau_{1i}, \quad \beta_{2i} = \beta_2 + \tau_{2i}, \quad \beta_{3i} = \beta_3 + \tau_{3i}, \quad \beta_{4i} = \beta_4 + \tau_{4i}, \quad \beta_{5i} = \beta_5 + \tau_{5i}, 
\]
\[
e_{ij} \sim N(0, \sigma_2^2), \quad \bm{\tau}_i \sim N(0, B), \quad i = 1, 2, \dots, N; \, j = 1, 2, \dots, n_i.
\]


- **GLMM:**

\[
\log\left(\dfrac{P(y_{ij}=1)}{P(y_{ij}=0)}\right) = (\alpha_0 + a_{0i}) + (\alpha_1 + a_{1i}) \cdot t_{ij} + (\alpha_2 + a_{2i}) \cdot x_{ij},
\]
where:
\[
\mathbf{a}_i \sim N(0, A), \quad x_{ij} \sim N(0, 1).
\]


## Installation

This project utilizes Monolix for some population model simulations. To use Monolix through R, ensure that `lixoftConnector` is properly set up on your machine, and that Monolix is installed. For more details on installation, refer to the [Lixoft documentation](https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/)

## Simulation Functions

The key functions in this repository simulate data based on the models mentioned above. They allow users to assess the performance of GLMM and NLME models under different conditions, including varying the correlation structure of random effects, the model parameters, and sample sizes, etc.

### 1. `SimFunDecay()`

This function simulates data for longitudinal outcomes using the NLME model (Wu and Ding, 1999):
\[
y_{ij} = \log_{10} \left( e^{P_{1i} - \lambda_{1i} t_{ij}} + e^{P_{2i}} \right) + e_{ij}, 
\]
where:
\[
P_{1i} = P_1 + b_{1i}, \quad P_{2i} = P_2 + b_{2i}, \quad \lambda_{1i} = \lambda_1 + b_{3i}, 
\]
\[
e_{ij} \sim N(0, \sigma_1^2), \quad \mathbf{b}_i \sim N(0, D).
\]

#### Arguments:
p1 = 17, lambda1 = 4, p2 = 3, D = diag(c(2, 0.1, 0.1)),
            sigma = 0.1, time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 20
- `p1`: the population parameter $P_1$.
- `lambda1`: the population parameter $\lambda_1$.
- `p2`: the population parameter $P_2$.
- `D`: the variance covariance matrix of the random effects $\mathbf{b}_i$. Any of the diagonal elements $D_{11}$, $D_{22}$, and $D_{33}$ can be zero, which means that the corresponding population parameter does not have a random effect.


#### Example Usage:
```r
data <- simulate_nlme(
  n_subjects = 100,
  n_obs = 10,
  theta = c(2, 0.5),
  sigma_u = 0.1,
  model = function(time, theta) {
    theta[1] * exp(-theta[2] * time)  # Exponential decay model
  }
)
```

### 1. `simulate_glmm()`

This function simulates data for binary outcomes using a GLMM with random intercepts and/or slopes.

#### Arguments:
- `n_subjects`: Number of subjects in the simulation.
- `n_obs`: Number of observations per subject.
- `beta`: Vector of fixed effect coefficients.
- `sigma_u`: Standard deviation of the random intercept.
- `rho`: Correlation between random intercept and slope.
- `family`: Family of the GLMM model. Default is `binomial`.

#### Example Usage:
```r
data <- simulate_glmm(
  n_subjects = 100,
  n_obs = 10,
  beta = c(0.5, -0.2),
  sigma_u = 1.0,
  rho = 0.2,
  family = binomial()
)
```


