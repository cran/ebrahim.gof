## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
#  # Install from GitHub
#  devtools::install_github("ebrahimkhaled/ebrahim.gof")
#  
#  # Load the package
#  library(ebrahim.gof)

## ----message=FALSE------------------------------------------------------------
library(ebrahim.gof)

## -----------------------------------------------------------------------------
# Simulate binary data
set.seed(123)
n <- 500
x <- rnorm(n)
linpred <- 0.5 + 1.2 * x
prob <- plogis(linpred)  # Convert to probabilities
y <- rbinom(n, 1, prob)

# Fit logistic regression
model <- glm(y ~ x, family = binomial())
predicted_probs <- fitted(model)

# Perform Ebrahim-Farrington test
result <- ef.gof(y, predicted_probs, G = 10)
print(result)

## -----------------------------------------------------------------------------
# Test with different numbers of groups
group_sizes <- c(4, 8, 10, 15, 20)
results <- data.frame(
  Groups = group_sizes,
  P_value = sapply(group_sizes, function(g) {
    ef.gof(y, predicted_probs, G = g)$p_value
  })
)
print(results)

## -----------------------------------------------------------------------------
# Hosmer-Lemeshow test (requires ResourceSelection package)
if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  library(ResourceSelection)
  
  # Perform both tests
  ef_result <- ef.gof(y, predicted_probs, G = 10)
  hl_result <- hoslem.test(y, predicted_probs, g = 10)
  
  # Compare results
  comparison <- data.frame(
    Test = c("Ebrahim-Farrington", "Hosmer-Lemeshow"),
    P_value = c(ef_result$p_value, hl_result$p.value),
    Test_Statistic = c(ef_result$Test_Statistic, hl_result$statistic)
  )
  print(comparison)
} else {
  cat("ResourceSelection package not available for comparison\n")
}

## -----------------------------------------------------------------------------
# Directed Ebrahim-Farrington test (takes the fitted model)
def.gof(model)                       # default poly3 basis
def.gof(model, basis = "ensemble")   # combine all three bases (Cauchy)

# Ensemble of the three DEF bases
def.ensemble.gof(model)
def.ensemble.gof(model, add_ef = TRUE)   # add the omnibus EF

## -----------------------------------------------------------------------------
run.all.gof(model, include_slow = FALSE)

## -----------------------------------------------------------------------------
# Function to simulate power under model misspecification
simulate_power <- function(n, beta_quad = 0.1, n_sims = 100, G = 10) {
  rejections_ef <- 0
  rejections_hl <- 0
  
  for (i in 1:n_sims) {
    # Generate data with quadratic term (true model)
    x <- runif(n, -2, 2)
    linpred_true <- 0 + x + beta_quad * x^2
    prob_true <- plogis(linpred_true)
    y <- rbinom(n, 1, prob_true)
    
    # Fit misspecified linear model (omitting quadratic term)
    model_mis <- glm(y ~ x, family = binomial())
    pred_probs <- fitted(model_mis)
    
    # Ebrahim-Farrington test
    ef_test <- ef.gof(y, pred_probs, G = G)
    if (ef_test$p_value < 0.05) rejections_ef <- rejections_ef + 1
    
    # Hosmer-Lemeshow test (if available)
    if (requireNamespace("ResourceSelection", quietly = TRUE)) {
      hl_test <- ResourceSelection::hoslem.test(y, pred_probs, g = G)
      if (hl_test$p.value < 0.05) rejections_hl <- rejections_hl + 1
    }
  }
  
  power_ef <- rejections_ef / n_sims
  power_hl <- if (requireNamespace("ResourceSelection", quietly = TRUE)) {
    rejections_hl / n_sims
  } else {
    NA
  }
  
  return(list(power_ef = power_ef, power_hl = power_hl))
}

# Calculate power for different sample sizes
sample_sizes <- c(200, 500, 1000)
power_results <- data.frame(
  n = sample_sizes,
  EbrahimFarrington_Power = sapply(sample_sizes, function(n) {
    simulate_power(n, beta_quad = 0.15, n_sims = 50)$power_ef
  })
)

if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  power_results$HosmerLemeshow_Power <- sapply(sample_sizes, function(n) {
    simulate_power(n, beta_quad = 0.15, n_sims = 50)$power_hl
  })
}

print(power_results)

## -----------------------------------------------------------------------------
# Simulate grouped data
set.seed(456)
n_groups <- 30
m_trials <- sample(5:20, n_groups, replace = TRUE)
x_grouped <- rnorm(n_groups)
prob_grouped <- plogis(0.2 + 0.8 * x_grouped)
y_grouped <- rbinom(n_groups, m_trials, prob_grouped)

# Create data frame and fit model
data_grouped <- data.frame(
  successes = y_grouped,
  trials = m_trials,
  x = x_grouped
)

model_grouped <- glm(
  cbind(successes, trials - successes) ~ x,
  data = data_grouped,
  family = binomial()
)

predicted_probs_grouped <- fitted(model_grouped)

# Original Farrington test for grouped data
result_grouped <- ef.gof(
  y_grouped,
  predicted_probs_grouped,
  model = model_grouped,
  m = m_trials,
  G = NULL  # No automatic grouping for original test
)

print(result_grouped)

