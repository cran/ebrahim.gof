## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
library(ebrahim.gof)

## ----data---------------------------------------------------------------------
data("gof_demo", package = "ebrahim.gof")
str(gof_demo)

fit <- glm(outcome ~ age + bmi + sex + treatment,
           data = gof_demo, family = binomial)

## ----battery------------------------------------------------------------------
set.seed(1)
battery <- run.all.gof(fit, include_slow = FALSE, install = "no")
battery

## ----grouped------------------------------------------------------------------
data("gof_demo_grouped", package = "ebrahim.gof")
fit_g <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo_grouped, family = binomial)
set.seed(1)
run.all.gof(fit_g, include_slow = FALSE, install = "no")

## ----directed-----------------------------------------------------------------
edge.gof(fit)                 # directed EDGE test (poly3 basis)
def.gof(fit, basis = "poly2") # a lower-order directed basis
ef.gof(fit)                   # the omnibus EF test

## ----ensemble-----------------------------------------------------------------
def.ensemble.gof(fit)

## -----------------------------------------------------------------------------
set.seed(4)
n <- 300; p <- 8
X <- matrix(rnorm(n * p), n, p)
y <- rbinom(n, 1, plogis(0.2 + X %*% c(0.9, -0.6, 0.4, rep(0, p - 3))))

# closed form: one fit, no resampling, and no Monte Carlo error in the p-value
calm.gof(X, y, lambda = 40)

