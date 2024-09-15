# Comparison of Computationally Efficient Approximate Methods for Nonlinear and Generalized Linear Mixed Effects Models

Generalized linear mixed models (GLMMs) and nonlinear mixed effects (NLME) models are popular in the analysis of longitudinal or clustered data. Statistical inference is typically based on likelihood methods. When the number of random effects in the models is large, the observed-data likelihood function involves high-dimensional and intractable integration, as these models are nonlinear in the (unobserved) random effects. "Exact" methods, such as Monte Carlo EM (MCEM) algorithms and numerical integration methods, can be computationally very intensive and may offer convergence issues. Computationally more efficient approximate methods, such as the stochastic approximation EM (SAEM) algorithm, linearization methods, or Laplace approximation methods, are therefore commonly used in practice. This repository contains the simulation code and functions used for evaluating the performance of these three methods for NLME and GLMM models. The SAEM, linearization, and Laplace approximation methods were implemented using the R packages `lixoftConnectors`, `nlme`, and `lme4`, respectively.

## Installation

This project utilizes Monolix for some population model simulations. To use Monolix through R, ensure that `lixoftConnector` is properly set up on your machine, and that Monolix is installed. For more details on installation, refer to the [Lixoft documentation](https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/).

## Simulation Functions

The key functions in this repository enable users to evaluate the performance of three approximation methods under various conditions, including different models, random effect correlation structures, model parameters, sample sizes, etc.

### 1. `SimFunDecay()`

This function simulates data for longitudinal outcomes using the NLME model (Wu and Ding, 1999):

$$y_{ij} = \log_{10} \left( e^{P_{1i} - \lambda_{1i} t_{ij}} + e^{P_{2i}} \right) + e_{ij}, $$

where:

$$P_{1i} = P_1 + b_{1i}, \quad P_{2i} = P_2 + b_{2i}, \quad \lambda_{1i} = \lambda_1 + b_{3i}, \quad e_{ij} \sim N(0, \sigma_1^2), \quad \mathbf{b}_i \sim N(0, D).$$

#### Arguments:

- `p1_t`: the population parameter $P_1$.
- `lambda1_t`: the population parameter $\lambda_1$.
- `p2_t`: the population parameter $P_2$.
- `D`: the variance-covariance matrix of the random effects. Any of the diagonal elements $D_{11}$, $D_{22}$, and $D_{33}$ can be zero, indicating that the corresponding population parameter does not have a random effect.
- `sigma`: the residual SD $\sigma_1$.
- `time`: a vector of time points.
- `num_sim`: the number of simulation iterations.
- `num_patient`: the number of simulated individuals.


#### Example Usage:

```r
SimFunDecay(p1_t = 17, lambda1_t = 4, p2_t = 3, D = diag(c(2, 0.1, 0.1)),
            sigma = 0.1, time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 20)
```

The result table will be stored in a folder named `Decay_`, followed by the parameter information.

### 2. `SimFunRebound()`

This function simulates data for longitudinal outcomes using the NLME model (Wang et al., 2020):

$$y_{ij} = \dfrac{\beta_{1i} t_{ij}}{t  + \exp{(\beta_{2i} - \beta_{3i} t_{ij})}} + \dfrac{\beta_{4i}}{1 + \exp( \beta_{5i} t_{ij})} + e_{ij}, $$

where:

$$\beta_{1i} = \beta_1 + \tau_{1i}, \quad \beta_{2i} = \beta_2 + \tau_{2i}, \quad \beta_{3i} = \beta_3 + \tau_{3i}, \quad \beta_{4i} = \beta_4 + \tau_{4i}, \quad \beta_{5i} = \beta_5 + \tau_{5i}, \quad e_{ij} \sim N(0, \sigma_2^2), \quad {\mathbf\tau}_i \sim N(0, B)$$


#### Arguments:

- `beta1_t`, `beta2_t`, `beta3_t`, `beta4_t`, and `beta5_t`: the population parameters $\beta_1$, $\beta_2$, $\beta_3$, $\beta_4$, and $\beta_5$.
- `B`: the variance-covariance matrix of the random effects. Any of the diagonal elements $B_{11}$, $B_{22}$, $B_{33}$, and $B_{44}$ can be zero, indicating that the corresponding population parameter does not have a random effect.
- `sigma`: the residual SD $\sigma_2$.
- `time`: a vector of time points.
- `num_sim`: the number of simulation iterations.
- `num_patient`: the number of simulated individuals.

#### Example Usage:

```r
SimFunRebound(beta1_t = 4, beta2_t = 5.7, beta3_t = 2.1, beta4_t = 1.9, beta5_t = 0.4,
              B = diag(c(0.2, 0.1, 0.1, 0, 0)), sigma = 0.1, 
              time = c(0.5, 2.9, 4.8, 7, 10.2), num_sim = 500, num_patient = 50)

```

The result table will be stored in a folder named `Rebound_`, followed by the parameter information.

### 3. `SimFunGLMM()`

This function simulates data for binary longitudinal outcomes using the GLMM:

$$\log\left(\dfrac{P(y_{ij}=1)}{P(y_{ij}=0)}\right) = (\alpha_0 + a_{0i}) + (\alpha_1 + a_{1i}) \cdot t_{ij} + (\alpha_2 + a_{2i}) \cdot x_{ij},$$

where:

$$ a_i \sim N(0, A), \quad x_{ij} \sim N(0, 1).$$


#### Arguments:

- `alpha0_t`, `alpha1_t`, and `alpha2_t`: the population parameter $\alpha_0$, $\alpha_1$, and $\alpha_2$.
- `A`: the variance covariance matrix of the random effects. Any of the diagonal elements $A_{11}$, $A_{22}$, and $A_{33}$ can be zero, indicating that the corresponding population parameter does not have a random effect.
- `time`: a vector of time points.
- `num_sim`: the number of simulation iterations.
- `num_patient`: the number of simulated individuals.


#### Example Usage:

```r
SimFunGLMM(alpha0_t = 8, alpha1_t = -2, alpha2_t = 2, A = diag(c(2, 0.5, 0)),
           time = c(0.5, 1.67, 3, 4.6, 6), num_sim = 500, num_patient = 50)
```

The result table will be stored in a folder named `GLMM_`, followed by the parameter information.



