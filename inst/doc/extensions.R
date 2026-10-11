## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 3.6)

## ----demo-fit-----------------------------------------------------------------
library("ebrahim.gof")
data("gof_demo", package = "ebrahim.gof")
fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo,
  family = binomial)

## ----grouped------------------------------------------------------------------
data("gof_demo_grouped", package = "ebrahim.gof")
fit_g <- glm(outcome ~ age + bmi + sex + treatment,
  data = gof_demo_grouped, family = binomial)
run.all.gof(fit_g, tests = c("Pearson", "Pearson (grouped)",
  "Deviance", "Deviance (grouped)", "McCullagh",
  "McCullagh (grouped)", "Osius-Rojek"))

## ----ef-----------------------------------------------------------------------
ef <- ef.gof(fit, G = "auto", reference = "normal")
ef
g <- attr(ef, "groups")
head(g[order(abs(g$c), decreasing = TRUE),
  c("group", "n", "observed", "expected", "pbar", "c")], 4)

## ----dealias------------------------------------------------------------------
loc_d <- localize.gof(fit, B = 199, seed = 1, dealias = TRUE)
loc_d

## ----stream-------------------------------------------------------------------
set.seed(2)
p_dev <- plogis(rnorm(5000, -1.5, 1))
s <- edge.stream(p_ref = p_dev, G = 10)
pv_day <- numeric(60)
for (day in 1:60) {
  lp <- rnorm(200, -1.5, 1)
  yy <- rbinom(200, 1, plogis(if (day <= 30) lp else 0.7 * lp))
  s <- update(s, yy, plogis(lp))
  pv_day[day] <- summary(s)$p_value
}
s
round(pv_day[c(30, 31, 32, 33, 35, 40)], 4)

## ----stream-plot, fig.height = 3.2, fig.alt = "Daily minus log10 p-value of edge.stream; flat before day 31, rising steeply after it."----
lp10 <- -log10(pmax(pv_day, 1e-300))
op <- par(mar = c(3.6, 3.9, 0.6, 0.6), mgp = c(2.4, 0.6, 0), las = 1)
plot(1:60, lp10, type = "n", xlab = "Day of deployment",
  ylab = "-log10 of the p-value to date", bty = "l")
rect(30.5, -1, 61, max(lp10) * 2, col = "#F2F2F2", border = NA)
abline(v = 30.5, lty = 2, col = "#D55E00")
abline(h = -log10(0.05), lty = 3, col = "grey40")
lines(1:60, lp10, type = "s", lwd = 1.6, col = "#0072B2")
par(op)

## ----ridge-data---------------------------------------------------------------
set.seed(4)
n <- 400
k <- 80
Xr <- matrix(rnorm(n * k), n, k)
yr <- rbinom(n, 1, plogis(-0.5 + Xr[, 1:5] %*% rep(0.5, 5)))
sum(yr)

## ----shrink-------------------------------------------------------------------
shrink.gof(Xr, yr, lambda = 40, B = 199, seed = 1,
  uncorrected = TRUE)

## ----calm---------------------------------------------------------------------
calm.gof(Xr, yr, lambda = 40, basis = "adaptive")

## ----lambda-scale-------------------------------------------------------------
a1 <- calm.gof(Xr, yr, lambda = 40, basis = "adaptive")
a2 <- calm.gof(Xr, yr, lambda = 0.1, lambda_scale = "glmnet",
  basis = "adaptive")
all.equal(a1, a2)

