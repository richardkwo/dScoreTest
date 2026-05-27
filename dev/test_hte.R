# Temporary test script for dst.hte
# Run from the dst package root: source("test_hte.R")

set.seed(NULL)

# ---- Packages ---------------------------------------------------------------
library(grf)

# ---- Source package files ---------------------------------------------------
source("R/dst.R")
source("R/hunt.R")
source("R/gam.R")
source("R/hte.R")

# ---- Data generation --------------------------------------------------------
# PLM: Y = g(Z) + D*f(Z) + eps
# Null (het.vars = NULL): f(Z) = beta (constant)
# Alternative:            f(Z) = beta + tau/sqrt(n) * sin(pi*Z[,1])
#
# Treatment: binary, P(D=1|Z) = sigmoid(Z[,1]/8 + 0.25*sin(pi*Z[,2]))
# Nuisance:  g(Z) = sigmoid((Z[,2]+Z[,3])/2) + Z[,1]

sigmoid <- function(x) 1 / (1 + exp(-x))

gen_hte_data <- function(n = 2000, p = 7, tau = 0) {
  Z <- matrix(runif(n * p, -1, 1), nrow = n,
              dimnames = list(NULL, paste0("Z", 1:p)))
  D <- rbinom(n, 1, sigmoid(Z[, 1] / 8 + 0.25 * sin(pi * Z[, 2])))
  g <- sigmoid((Z[, 2] + Z[, 3]) / 2) + Z[, 1]
  f <- 0.75 + (tau / sqrt(n)) * sin(pi * Z[, 1])
  y <- g + D * f + sqrt(0.25 * (1 + Z[, 2]^2)) * rnorm(n)
  list(y = y, Z = Z, D = D)
}

# ---- Helpers ----------------------------------------------------------------
pval_label <- function(p) {
  if (p < 0.01) "** (p < 0.01)"
  else if (p < 0.05) "*  (p < 0.05)"
  else sprintf("   (p = %.3f, not significant)", p)
}

cat_result <- function(label, p) {
  cat(sprintf("%-60s p = %.4f  %s\n", label, p, pval_label(p)))
}

# =============================================================================
cat("\n====  dst.hte  (null: constant effect)  ====\n\n")

d0 <- gen_hte_data(n = 2000, p = 7, tau = 0)
d5 <- gen_hte_data(n = 2000, p = 7, tau = 5)

p1 <- dst.hte(d0$y, d0$Z, d0$D)
cat_result("null  (tau=0)", p1)

p2 <- dst.hte(d5$y, d5$Z, d5$D)
cat_result("alt   (tau=5)", p2)

# =============================================================================
cat("\n====  dst.hte  (hunt styles)  ====\n\n")

p3 <- dst.hte(d5$y, d5$Z, d5$D, hunt.style = "WLS")
cat_result("alt (tau=5, WLS)", p3)

p4 <- dst.hte(d5$y, d5$Z, d5$D, hunt.style = "vanilla")
cat_result("alt (tau=5, vanilla)", p4)

# =============================================================================
cat("\n====  dst.hte  (het.vars: test specific variables)  ====\n\n")
cat("True heterogeneity in Z1. het.vars=1 tests Z1 specifically.\n\n")

# het.vars = 1: null allows heterogeneity in Z2..Z7, test asks about Z1
p5 <- dst.hte(d0$y, d0$Z, d0$D, het.vars = 1)
cat_result("null  (tau=0, het.vars=1)", p5)

p6 <- dst.hte(d5$y, d5$Z, d5$D, het.vars = 1)
cat_result("alt   (tau=5, het.vars=1)", p6)

# het.vars = 2: null allows heterogeneity in Z1,Z3..Z7, test asks about Z2
# Z2 has no true heterogeneity — should not reject
p7 <- dst.hte(d5$y, d5$Z, d5$D, het.vars = 2)
cat_result("alt   (tau=5, het.vars=2, no true het in Z2)", p7)

# =============================================================================
cat("\n====  dst.hte  (return.hunted)  ====\n\n")

res <- dst.hte(d5$y, d5$Z, d5$D, return.hunted = TRUE)
cat(sprintf("p-value: %.4f\n", res$p.value))
cat(sprintf("h.hat is a function: %s\n", is.function(res$h.hat)))

# =============================================================================
cat("\n====  dst.hte  (multi.split)  ====\n\n")

p8 <- dst.hte(d0$y, d0$Z, d0$D, multi.split = TRUE, B = 5)
cat_result("null  (tau=0, multi.split)", p8)

p9 <- dst.hte(d5$y, d5$Z, d5$D, multi.split = TRUE, B = 5)
cat_result("alt   (tau=5, multi.split)", p9)

# =============================================================================
cat("\n====  Variable importance via het.vars  ====\n\n")
cat("DGP: CATE = 0.75 + tau/sqrt(n) * (Z1 + Z2 + Z3) / sqrt(3)\n")
cat("Each variable contributes equally -- importance scores should be roughly equal.\n\n")

# 3-variable DGP with equal contributions to CATE
n.vi  <- 2000
tau.vi <- 8
Z.vi <- matrix(runif(n.vi * 3, -1, 1), nrow = n.vi,
               dimnames = list(NULL, paste0("Z", 1:3)))
D.vi <- rbinom(n.vi, 1, sigmoid(Z.vi[, 1] / 4))
f.vi <- 0.75 + (tau.vi / sqrt(n.vi)) * (Z.vi[, 1] + Z.vi[, 2] + Z.vi[, 3]) / sqrt(3)
y.vi <- sigmoid(Z.vi[, 2]) + D.vi * f.vi + 0.5 * rnorm(n.vi)

# Run DST for each variable
importance <- data.frame(variable = paste0("Z", 1:3), p.value = NA_real_)
for (i in 1:3) {
  importance$p.value[i] <- dst.hte(y.vi, Z.vi, D.vi, het.vars = i)
}

importance$neg.log.p <- -log(importance$p.value)
importance <- importance[order(importance$p.value), ]

cat("Variable importance (ranked by p-value):\n\n")
cat(sprintf("  %-10s  %-10s  %-12s\n", "Variable", "p-value", "-log(p-value)"))
cat(sprintf("  %-10s  %-10s  %-12s\n", "--------", "-------", "------------"))
for (i in seq_len(nrow(importance))) {
  cat(sprintf("  %-10s  %-10.4f  %-12.3f\n",
              importance$variable[i],
              importance$p.value[i],
              importance$neg.log.p[i]))
}
cat("\nAll three variables should have similar importance.\n")

cat("\n---- Done ----\n")
