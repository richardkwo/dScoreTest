# Temporary test script for dst.qgam
# Run from the dst package root: source("test_qgam.R")

set.seed(NULL)

# ---- Packages ---------------------------------------------------------------
library(mgcv)
library(qgam)
library(grf)
library(MASS)

# ---- Source package files ---------------------------------------------------
source("R/dst.R")
source("R/hunt.R")
source("R/gam.R")   # needed for .gam_fit used inside qgam.R
source("R/qgam.R")

# ---- Data generation --------------------------------------------------------
# Null:        y = sin(X1) + scale(X) * eps
# Alternative: y = sin(X1) + tau/sqrt(n) * X1*X3 + scale(X) * eps
# scale(X) = 0.25 * (1 + X2^2)  -- heteroskedastic
# X independent uniform on [-1, 1]

gen_qgam_data <- function(n = 1000, p = 5, tau = 0, err = "t") {
  X <- matrix(runif(n * p, -1, 1), nrow = n,
              dimnames = list(NULL, paste0("X", 1:p)))
  scale <- 0.25 * (1 + X[, 2]^2)
  eps   <- switch(err,
                  t         = rt(n, df = 3),
                  exp       = rexp(n) - 1,   # centred
                  gaussian  = rnorm(n))
  y <- sin(X[, 1]) + (tau / sqrt(n)) * X[, 1] * X[, 3] + scale * eps
  list(y = y, X = X)
}

# ---- Helpers ----------------------------------------------------------------
pval_label <- function(p) {
  if (p < 0.01) "** (p < 0.01)"
  else if (p < 0.05) "*  (p < 0.05)"
  else sprintf("   (p = %.3f, not significant)", p)
}

cat_result <- function(label, p) {
  cat(sprintf("%-55s p = %.4f  %s\n", label, p, pval_label(p)))
}

# =============================================================================
cat("\n====  dst.qgam  (qu = 0.5, t errors)  ====\n\n")

d0 <- gen_qgam_data(n = 1000, p = 5, tau = 0,  err = "t")
d5 <- gen_qgam_data(n = 1000, p = 5, tau = 10, err = "t")

p1 <- dst.qgam(d0$y, d0$X, qu = 0.5)
cat_result("null  (tau=0)",  p1)

p2 <- dst.qgam(d5$y, d5$X, qu = 0.5)
cat_result("alt   (tau=10)", p2)

# =============================================================================
cat("\n====  dst.qgam  (qu = 0.7, t errors)  ====\n\n")

p3 <- dst.qgam(d0$y, d0$X, qu = 0.7)
cat_result("null  (tau=0)",  p3)

p4 <- dst.qgam(d5$y, d5$X, qu = 0.7)
cat_result("alt   (tau=10)", p4)

# =============================================================================
cat("\n====  dst.qgam  (qu = 0.7, exponential errors)  ====\n\n")

de0 <- gen_qgam_data(n = 1000, p = 5, tau = 0,  err = "exp")
de5 <- gen_qgam_data(n = 1000, p = 5, tau = 10, err = "exp")

p5 <- dst.qgam(de0$y, de0$X, qu = 0.7)
cat_result("null  (tau=0)",  p5)

p6 <- dst.qgam(de5$y, de5$X, qu = 0.7)
cat_result("alt   (tau=10)", p6)

# =============================================================================
cat("\n====  dst.qgam  (return.hunted)  ====\n\n")

res <- dst.qgam(d5$y, d5$X, qu = 0.7, return.hunted = TRUE)
cat(sprintf("p-value: %.4f\n", res$p.value))
cat(sprintf("h.hat is a function: %s\n", is.function(res$h.hat)))
cat(sprintf("h.hat on first 3 obs: %s\n",
            paste(round(res$h.hat(d5$X[1:3, ]), 3), collapse = ", ")))

cat("\n---- Done ----\n")
