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

## ----battery-summary----------------------------------------------------------
summary(battery)

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

## ----edge-plot, fig.height = 3.6----------------------------------------------
e <- edge.gof(fit)
plot(e)

## ----external-----------------------------------------------------------------
dev <- gof_demo[1:400, ]
val <- gof_demo[401:800, ]
covs <- c("age", "bmi", "sex", "treatment")
fit_dev <- glm(outcome ~ age + bmi + sex + treatment, data = dev, family = binomial)
p_val <- predict(fit_dev, newdata = val, type = "response")   # frozen
run.all.external(val$outcome, p_val)

## ----external-edge, fig.height = 3.6------------------------------------------
e_ext <- edge.gof(y = val$outcome, predicted_probs = p_val, external = TRUE)
e_ext
plot(e_ext, scale = "rate")

## ----localize, fig.height = 3.4-----------------------------------------------
loc <- localize.gof(fit, B = 199, seed = 1)
loc
plot(loc)

## ----localize-external, fig.height = 3.4--------------------------------------
loc_ext <- localize.external(val$outcome, p_val, val[, covs], M = 199, seed = 1)
plot(loc_ext)

## ----stream-------------------------------------------------------------------
set.seed(2)
s <- edge.stream(p_ref = p_val, G = 10)
for (day in 1:60) {                       # 100 patients a day
  lp <- rnorm(100, -1, 1)
  shift <- if (day > 30) 0.5 else 0      # from day 31 the event rate drifts upwards
  s <- update(s, rbinom(100, 1, plogis(lp + shift)), plogis(lp))
  if (day %in% c(30, 60)) print(s)
}

## -----------------------------------------------------------------------------
set.seed(4)
n <- 300; p <- 8
X <- matrix(rnorm(n * p), n, p)
y <- rbinom(n, 1, plogis(0.2 + X %*% c(0.9, -0.6, 0.4, rep(0, p - 3))))

# closed form: one fit, no resampling, and no Monte Carlo error in the p-value
calm.gof(X, y, lambda = 40)

