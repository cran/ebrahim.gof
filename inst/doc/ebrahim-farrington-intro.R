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
set.seed(2026)

power_paired <- function(n, beta_quad = 0.15, n_sims = 200, G = 10) {
  have_rs <- requireNamespace("ResourceSelection", quietly = TRUE)
  rej <- c(EF = 0, HL = 0)
  for (i in seq_len(n_sims)) {
    x    <- runif(n, -2, 2)
    y    <- rbinom(n, 1, plogis(x + beta_quad * x^2))
    fit  <- glm(y ~ x, family = binomial())      # misspecified: no quadratic term
    p    <- fitted(fit)

    # both tests see this same dataset
    if (ef.gof(y, p, G = G)$p_value < 0.05) rej["EF"] <- rej["EF"] + 1
    if (have_rs &&
        ResourceSelection::hoslem.test(y, p, g = G)$p.value < 0.05) rej["HL"] <- rej["HL"] + 1
  }
  c(n = n, EF = unname(rej["EF"]) / n_sims,
    HL = if (have_rs) unname(rej["HL"]) / n_sims else NA_real_)
}

power_results <- as.data.frame(do.call(rbind, lapply(c(200, 500, 1000), power_paired)))
print(power_results, row.names = FALSE)

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

