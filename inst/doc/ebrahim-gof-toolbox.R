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

